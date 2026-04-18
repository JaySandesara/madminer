#!/usr/bin/env python
"""
Train a parameterized likelihood-ratio estimator on combined_samples.h5.

The network learns  log r(x | theta, nu) = log p(x | theta, nu) / p(x | theta_ref, nu_ref),
with theta_ref = (0, 0) and nu_ref = (0, 0).

Training targets come for free from the morphing + nuisance model already
stored per event:
    w(x | theta, nu) = (morphing_matrix @ t(theta)) . w_benchmarks(x)
                       * exp( a . nu + b . nu^2 )

For each minibatch we draw random (theta, nu) from priors and regress the NN
output onto the truth log-ratio (ALICES-style direct regression, no classifier).

Usage:
    python train_parameterized_carl.py \\
        --input /staging/jsandesara/madminer/semi_parametric/combined_samples.h5 \\
        --outdir ./carl_model \\
        --epochs 20

Requires: numpy, h5py, torch, matplotlib.
"""

import argparse
import os
import time
from pathlib import Path

import h5py
import numpy as np
import torch
import torch.nn as nn


# ======================================================================
# Data loading and weight computation
# ======================================================================

def load_combined(path):
    with h5py.File(path, "r") as f:
        data = {
            "x": f["observables"][()],
            "w_bm": f["weights_benchmarks"][()],
            "nuis_a": f["nuisance_a"][()],        # (n_nuis, n_events)
            "nuis_b": f["nuisance_b"][()],
            "morph_mat": f["morphing_matrix"][()],       # (n_bm, n_comp)
            "morph_comp": f["morphing_components"][()],  # (n_comp, n_params)
            "bm_vals": f["benchmark_values"][()],
            "param_names": [s.decode() for s in f["parameter_names"][()]],
            "obs_names":   [s.decode() for s in f["observable_names"][()]],
            "nu_names":    [s.decode() for s in f["nuisance_parameter_names"][()]],
        }
    return data


def t_vector_batch(theta, comp):
    """Polynomial basis t_j(theta) for a batch of theta.

    theta:  (B, n_params)
    comp:   (n_comp, n_params)  — integer powers per monomial
    Returns (B, n_comp).
    """
    # theta[:, None, :]**comp[None, :, :] -> (B, n_comp, n_params) -> prod axis=-1
    return np.prod(theta[:, None, :] ** comp[None, :, :], axis=-1)


def weights_at_theta_nu_batch(theta, nu, x_indices, data):
    """Compute per-event weight w(x_i | theta_i, nu_i) for each sample in a batch.

    theta: (B, n_params)
    nu:    (B, n_nuis)
    x_indices: (B,) indices into the full event array
    Returns: (B,) per-sample weight (in pb, same units as weights_benchmarks).
    """
    # Morphing coefficients per benchmark for each theta in the batch
    t_vec = t_vector_batch(theta, data["morph_comp"])                 # (B, n_comp)
    morph_coeffs = t_vec @ data["morph_mat"].T                         # (B, n_bm)

    # Physics weight for the picked event at each batch element's theta
    w_bm_picked = data["w_bm"][x_indices]                              # (B, n_bm)
    w_phys = np.sum(w_bm_picked * morph_coeffs, axis=1)                # (B,)

    # Nuisance factor: exp( a . nu + b . nu^2 )
    a_picked = data["nuis_a"][:, x_indices].T                          # (B, n_nuis)
    b_picked = data["nuis_b"][:, x_indices].T                          # (B, n_nuis)
    log_factor = np.sum(a_picked * nu + b_picked * (nu**2), axis=1)    # (B,)

    return w_phys * np.exp(log_factor)


def weights_at_theta_nu_full(theta, nu, data):
    """Same as batch version but for ALL events at a single (theta, nu) point."""
    t_vec = np.prod(theta ** data["morph_comp"], axis=1)               # (n_comp,)
    morph_coeffs = data["morph_mat"] @ t_vec                           # (n_bm,)
    w_phys = data["w_bm"] @ morph_coeffs                               # (n_events,)
    log_factor = data["nuis_a"].T @ nu + data["nuis_b"].T @ (nu**2)    # (n_events,)
    return w_phys * np.exp(log_factor)


# ======================================================================
# Sampling
# ======================================================================

class PriorSampler:
    """Draws (theta, nu) from simple priors:
       theta ~ Uniform(-5, 5) per POI  (within generated benchmark grid)
       nu    ~ Gaussian(0, 1)          (clipped to [-3, 3])
    """

    def __init__(self, n_params, n_nuis, theta_range=(-5.0, 5.0), nu_sigma=1.0, nu_clip=3.0):
        self.n_params = n_params
        self.n_nuis = n_nuis
        self.theta_range = theta_range
        self.nu_sigma = nu_sigma
        self.nu_clip = nu_clip
        self.rng = np.random.default_rng()

    def sample(self, n):
        theta = self.rng.uniform(
            self.theta_range[0], self.theta_range[1], size=(n, self.n_params)
        )
        nu = np.clip(
            self.rng.normal(0.0, self.nu_sigma, size=(n, self.n_nuis)),
            -self.nu_clip, self.nu_clip,
        )
        return theta.astype(np.float32), nu.astype(np.float32)


def build_training_batch(batch_size, data, prior, rng):
    """Build (inputs, log_r_target) for one training minibatch.

    - Pick random events (uniform over all events)
    - Sample random (theta_num, nu_num) from prior
    - theta_den = 0, nu_den = 0 (fixed reference)
    - Compute log_r = log(w_num) - log(w_den) per event
    - Return (x, theta_num, nu_num, log_r) as float32 tensors.
    """
    n_events = data["x"].shape[0]
    idx = rng.integers(0, n_events, size=batch_size)
    theta_num, nu_num = prior.sample(batch_size)

    w_num = weights_at_theta_nu_batch(theta_num, nu_num, idx, data)
    theta_den = np.zeros_like(theta_num)
    nu_den = np.zeros_like(nu_num)
    w_den = weights_at_theta_nu_batch(theta_den, nu_den, idx, data)

    # Filter out problematic events (negative/zero weights from interference)
    good = (w_num > 0) & (w_den > 0) & np.isfinite(w_num) & np.isfinite(w_den)
    if not np.any(good):
        return None
    idx = idx[good]
    theta_num = theta_num[good]
    nu_num = nu_num[good]
    w_num = w_num[good]
    w_den = w_den[good]

    log_r = np.log(w_num) - np.log(w_den)
    # Clip extreme values to stabilize training
    log_r = np.clip(log_r, -15.0, 15.0)

    x = data["x"][idx].astype(np.float32)
    inp = np.concatenate([x, theta_num, nu_num], axis=1).astype(np.float32)
    return inp, log_r.astype(np.float32)


# ======================================================================
# Model
# ======================================================================

class ParamRatioMLP(nn.Module):
    def __init__(self, n_in, hidden=(256, 256, 256)):
        super().__init__()
        layers = []
        prev = n_in
        for h in hidden:
            layers += [nn.Linear(prev, h), nn.Tanh()]
            prev = h
        layers += [nn.Linear(prev, 1)]
        self.net = nn.Sequential(*layers)

    def forward(self, x):
        return self.net(x).squeeze(-1)


# ======================================================================
# Training loop
# ======================================================================

def standardize(data):
    """Return per-observable mean/std for input normalization."""
    mean = data["x"].mean(axis=0).astype(np.float32)
    std = data["x"].std(axis=0).astype(np.float32)
    std[std == 0] = 1.0
    return mean, std


def train(args):
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Device: {device}")

    print(f"Loading {args.input}")
    data = load_combined(args.input)
    n_events = data["x"].shape[0]
    n_obs = data["x"].shape[1]
    n_params = data["morph_comp"].shape[1]
    n_nuis = data["nuis_a"].shape[0]
    print(f"  events: {n_events:,}")
    print(f"  observables: {n_obs}  parameters: {n_params}  nuisances: {n_nuis}")
    print(f"  params: {data['param_names']}  nuis: {data['nu_names']}")

    # Input standardization
    x_mean, x_std = standardize(data)
    data["x"] = (data["x"] - x_mean) / x_std

    prior = PriorSampler(n_params, n_nuis)
    model = ParamRatioMLP(n_in=n_obs + n_params + n_nuis, hidden=tuple(args.hidden)).to(device)
    opt = torch.optim.Adam(model.parameters(), lr=args.lr)
    sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, T_max=args.epochs)
    loss_fn = nn.MSELoss()

    rng = np.random.default_rng(42)
    steps_per_epoch = max(1, n_events // args.batch_size)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    history = {"train_loss": [], "val_loss": []}

    val_fixed = []
    for _ in range(8):
        b = build_training_batch(args.batch_size, data, prior, rng)
        if b is not None:
            val_fixed.append(b)

    print(f"Training {args.epochs} epochs x {steps_per_epoch} steps/epoch, batch={args.batch_size}")

    for epoch in range(args.epochs):
        t0 = time.time()
        model.train()
        losses = []
        for step in range(steps_per_epoch):
            b = build_training_batch(args.batch_size, data, prior, rng)
            if b is None:
                continue
            inp, log_r = b
            inp_t = torch.from_numpy(inp).to(device)
            tgt_t = torch.from_numpy(log_r).to(device)
            pred = model(inp_t)
            loss = loss_fn(pred, tgt_t)
            opt.zero_grad()
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 5.0)
            opt.step()
            losses.append(loss.item())

        sched.step()
        train_loss = float(np.mean(losses))

        # Validation on fixed minibatches
        model.eval()
        vlosses = []
        with torch.no_grad():
            for inp, log_r in val_fixed:
                inp_t = torch.from_numpy(inp).to(device)
                tgt_t = torch.from_numpy(log_r).to(device)
                pred = model(inp_t)
                vlosses.append(loss_fn(pred, tgt_t).item())
        val_loss = float(np.mean(vlosses))

        history["train_loss"].append(train_loss)
        history["val_loss"].append(val_loss)
        dt = time.time() - t0
        print(f"  epoch {epoch+1:3d}/{args.epochs}  "
              f"train={train_loss:.4f}  val={val_loss:.4f}  "
              f"lr={sched.get_last_lr()[0]:.2e}  [{dt:.1f}s]")

    # Save everything
    torch.save(
        {
            "model_state_dict": model.state_dict(),
            "hidden": tuple(args.hidden),
            "n_in": n_obs + n_params + n_nuis,
            "n_obs": n_obs,
            "n_params": n_params,
            "n_nuis": n_nuis,
            "param_names": data["param_names"],
            "nu_names": data["nu_names"],
            "obs_names": data["obs_names"],
            "x_mean": x_mean,
            "x_std": x_std,
            "history": history,
        },
        outdir / "model.pt",
    )
    print(f"Saved model to {outdir / 'model.pt'}")

    # Plot
    try:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(8, 4))
        ax.plot(history["train_loss"], label="train")
        ax.plot(history["val_loss"], label="val")
        ax.set_xlabel("epoch")
        ax.set_ylabel("MSE on log r")
        ax.set_yscale("log")
        ax.legend()
        ax.set_title("Parameterized CARL training")
        plt.tight_layout()
        plt.savefig(outdir / "loss.png", dpi=130)
        print(f"Saved loss curve to {outdir / 'loss.png'}")
    except Exception as e:
        print(f"Skipped loss plot: {e}")

    # Closure check on a held-out random (theta, nu) scan
    print("\nClosure check (random 5k events, 1 random theta/nu per event):")
    model.eval()
    with torch.no_grad():
        n_test = 5000
        idx = rng.integers(0, n_events, size=n_test)
        theta_t, nu_t = prior.sample(n_test)
        w_num = weights_at_theta_nu_batch(theta_t, nu_t, idx, data)
        w_den = weights_at_theta_nu_batch(
            np.zeros_like(theta_t), np.zeros_like(nu_t), idx, data
        )
        good = (w_num > 0) & (w_den > 0)
        log_r_true = np.log(w_num[good]) - np.log(w_den[good])
        log_r_true = np.clip(log_r_true, -15, 15)
        x = data["x"][idx[good]].astype(np.float32)
        inp = np.concatenate([x, theta_t[good], nu_t[good]], axis=1).astype(np.float32)
        pred = model(torch.from_numpy(inp).to(device)).cpu().numpy()
        res = pred - log_r_true
        print(f"  true:  mean={log_r_true.mean():.3f}  std={log_r_true.std():.3f}")
        print(f"  pred:  mean={pred.mean():.3f}  std={pred.std():.3f}")
        print(f"  resid: mean={res.mean():.3f}  std={res.std():.3f}")


# ======================================================================
# Main
# ======================================================================

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--input",  default="/staging/jsandesara/madminer/semi_parametric/combined_samples.h5")
    p.add_argument("--outdir", default="./carl_model")
    p.add_argument("--epochs", type=int, default=20)
    p.add_argument("--batch-size", type=int, default=1024)
    p.add_argument("--lr", type=float, default=1e-3)
    p.add_argument("--hidden", type=int, nargs="+", default=[256, 256, 256])
    return p.parse_args()


if __name__ == "__main__":
    train(parse_args())
