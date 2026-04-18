"""
Aggregate multiple semi_parametric_samples_batch*.h5 files and optionally convert to ROOT.

Usage:
    python aggregate_and_convert.py --input-dir data/ --output data/combined_samples.h5
    python aggregate_and_convert.py --input-dir data/ --output data/combined_samples.h5 --root

Expects files matching: semi_parametric_samples_batch*.h5
Each file must have the same observables, benchmarks, and morphing setup.

Requires: h5py, numpy, uproot (only if --root)
"""

import argparse
import glob
import os
import sys

import h5py
import numpy as np


def load_batch(path):
    """Load a single batch file and return all arrays + metadata."""
    with h5py.File(path, "r") as f:
        data = {
            "observables": f["observables"][()],
            "weights_benchmarks": f["weights_benchmarks"][()],
            "nuisance_a": f["nuisance_a"][()],
            "nuisance_b": f["nuisance_b"][()],
            "sampling_ids": f["sampling_ids"][()],
        }

        # Load nuisance weight columns from group
        data["nuisance_weights"] = {}
        if "nuisance_weights" in f:
            for key in f["nuisance_weights"]:
                data["nuisance_weights"][key] = f["nuisance_weights"][key][()]

        # Load dense scale weights (per-event raw muF variation weights)
        data["scale_weights_dense"] = {}
        if "scale_weights_dense" in f:
            for syst_name in f["scale_weights_dense"]:
                data["scale_weights_dense"][syst_name] = {}
                for muf_key in f[f"scale_weights_dense/{syst_name}"]:
                    data["scale_weights_dense"][syst_name][muf_key] = \
                        f[f"scale_weights_dense/{syst_name}/{muf_key}"][()]

        meta = {
            "morphing_matrix": f["morphing_matrix"][()],
            "morphing_components": f["morphing_components"][()],
            "benchmark_names": f["benchmark_names"][()],
            "benchmark_values": f["benchmark_values"][()],
            "parameter_names": f["parameter_names"][()],
            "observable_names": f["observable_names"][()],
        }
        if "nuisance_parameter_names" in f:
            meta["nuisance_parameter_names"] = f["nuisance_parameter_names"][()]
        if "systematics_names" in f:
            meta["systematics_names"] = f["systematics_names"][()]
    return data, meta


def validate_compatible(meta_ref, meta_new, path_ref, path_new):
    """Check that two batch files have consistent metadata."""
    for key in ["benchmark_names", "parameter_names", "observable_names"]:
        if not np.array_equal(meta_ref[key], meta_new[key]):
            print(f"ERROR: {key} mismatch between {path_ref} and {path_new}")
            sys.exit(1)
    if not np.allclose(meta_ref["benchmark_values"], meta_new["benchmark_values"]):
        print(f"ERROR: benchmark_values mismatch between {path_ref} and {path_new}")
        sys.exit(1)
    if not np.allclose(meta_ref["morphing_matrix"], meta_new["morphing_matrix"]):
        print(f"ERROR: morphing_matrix mismatch between {path_ref} and {path_new}")
        sys.exit(1)


def aggregate(input_dir, output_path):
    """Find all batch files, validate, concatenate, and save."""
    pattern = os.path.join(input_dir, "semi_parametric_samples_batch*.h5")
    files = sorted(glob.glob(pattern))

    # Exclude the output file if it already exists in the same directory
    files = [f for f in files if os.path.abspath(f) != os.path.abspath(output_path)]

    if not files:
        print(f"No files matching {pattern}")
        sys.exit(1)

    print(f"Found {len(files)} batch files:")

    # Load first batch as reference
    all_data = []
    ref_data, ref_meta = load_batch(files[0])
    all_data.append(ref_data)
    print(f"  {files[0]}: {ref_data['observables'].shape[0]} events")

    # Load and validate remaining batches
    for path in files[1:]:
        data, meta = load_batch(path)
        validate_compatible(ref_meta, meta, files[0], path)
        all_data.append(data)
        print(f"  {path}: {data['observables'].shape[0]} events")

    # Concatenate event-level arrays
    combined = {}
    for key in ["observables", "weights_benchmarks", "sampling_ids"]:
        combined[key] = np.concatenate([d[key] for d in all_data], axis=0)

    # Nuisance coefficients: concat along event axis (axis=1)
    combined["nuisance_a"] = np.concatenate([d["nuisance_a"] for d in all_data], axis=1)
    combined["nuisance_b"] = np.concatenate([d["nuisance_b"] for d in all_data], axis=1)

    # Nuisance weight columns
    nuis_keys = sorted(all_data[0]["nuisance_weights"].keys())
    combined["nuisance_weights"] = {}
    for nk in nuis_keys:
        combined["nuisance_weights"][nk] = np.concatenate(
            [d["nuisance_weights"][nk] for d in all_data], axis=0
        )

    # Dense scale weights (aligned with events)
    combined["scale_weights_dense"] = {}
    if all_data[0]["scale_weights_dense"]:
        for syst_name, muf_dict in all_data[0]["scale_weights_dense"].items():
            combined["scale_weights_dense"][syst_name] = {}
            for muf_key in muf_dict:
                combined["scale_weights_dense"][syst_name][muf_key] = np.concatenate(
                    [d["scale_weights_dense"][syst_name][muf_key] for d in all_data], axis=0
                )

    n_total = combined["observables"].shape[0]
    print(f"\nTotal: {n_total} events")

    # Renormalize weights so the aggregate sum equals the per-batch mean sum.
    # Each batch's weights are scaled so that sum(weight_sm) ~ sigma(sm) (first physics
    # benchmark column). Concatenating N batches naively gives sum = N × sigma(sm);
    # we rescale back to ~sigma(sm).  nuisance_a/b are log-ratios and are left alone.
    per_batch_sum = np.array([d["weights_benchmarks"][:, 0].sum() for d in all_data])
    target_sum = per_batch_sum.mean()
    current_sum = combined["weights_benchmarks"][:, 0].sum()
    scale = target_sum / current_sum
    print(f"Renormalizing weights: target_sum={target_sum:.4e}, current_sum={current_sum:.4e}, scale={scale:.4e}")
    combined["weights_benchmarks"] = combined["weights_benchmarks"] * scale
    for nk in combined["nuisance_weights"]:
        combined["nuisance_weights"][nk] = combined["nuisance_weights"][nk] * scale
    for syst_name, muf_dict in combined["scale_weights_dense"].items():
        for muf_key in muf_dict:
            muf_dict[muf_key] = muf_dict[muf_key] * scale

    # Shuffle
    rng = np.random.default_rng(42)
    perm = rng.permutation(n_total)
    for key in combined:
        if key == "nuisance_weights":
            for nk in combined[key]:
                combined[key][nk] = combined[key][nk][perm]
        elif key == "scale_weights_dense":
            for syst_name, muf_dict in combined[key].items():
                for muf_key in muf_dict:
                    muf_dict[muf_key] = muf_dict[muf_key][perm]
        elif combined[key].ndim == 1:
            combined[key] = combined[key][perm]
        elif key in ("nuisance_a", "nuisance_b"):
            combined[key] = combined[key][:, perm]
        else:
            combined[key] = combined[key][perm]

    # Write combined HDF5
    with h5py.File(output_path, "w") as f:
        for key, arr in combined.items():
            if key == "nuisance_weights":
                grp = f.create_group("nuisance_weights")
                for nk, nv in arr.items():
                    grp.create_dataset(nk, data=nv)
            elif key == "scale_weights_dense":
                if arr:
                    grp = f.create_group("scale_weights_dense")
                    for syst_name, muf_dict in arr.items():
                        sgrp = grp.create_group(syst_name)
                        for muf_key, w in muf_dict.items():
                            sgrp.create_dataset(muf_key, data=w)
            else:
                f.create_dataset(key, data=arr)
        for key, arr in ref_meta.items():
            f.create_dataset(key, data=arr)

    print(f"Saved combined file to {output_path}")
    return output_path


def convert_to_root(h5_path, root_path):
    """Convert combined HDF5 to ROOT format."""
    try:
        import uproot
    except ImportError:
        print("uproot is required for ROOT conversion: pip install uproot")
        sys.exit(1)

    with h5py.File(h5_path, "r") as f:
        obs = f["observables"][()]
        w_bench = f["weights_benchmarks"][()]
        a_coeffs = f["nuisance_a"][()]
        b_coeffs = f["nuisance_b"][()]
        sampling_ids = f["sampling_ids"][()]
        morph_mat = f["morphing_matrix"][()]
        morph_comp = f["morphing_components"][()]
        bm_names = [s.decode() for s in f["benchmark_names"][()]]
        bm_vals = f["benchmark_values"][()]
        p_names = [s.decode() for s in f["parameter_names"][()]]
        o_names = [s.decode() for s in f["observable_names"][()]]

        # Nuisance weights
        nuis_weights = {}
        if "nuisance_weights" in f:
            for key in f["nuisance_weights"]:
                nuis_weights[key] = f["nuisance_weights"][key][()]

        nuis_param_names = []
        if "nuisance_parameter_names" in f:
            nuis_param_names = [s.decode() for s in f["nuisance_parameter_names"][()]]
        syst_names = []
        if "systematics_names" in f:
            syst_names = [s.decode() for s in f["systematics_names"][()]]

        # Dense scale weights (per-event raw muF variation)
        dense_weights = {}
        if "scale_weights_dense" in f:
            for syst_name in f["scale_weights_dense"]:
                dense_weights[syst_name] = {}
                for muf_key in f[f"scale_weights_dense/{syst_name}"]:
                    dense_weights[syst_name][muf_key] = \
                        f[f"scale_weights_dense/{syst_name}/{muf_key}"][()]

    # -- Events tree --
    event_data = {}

    # Observables
    for i, name in enumerate(o_names):
        event_data[name] = obs[:, i].astype(np.float64)

    # ME weights at each physical benchmark
    for i, name in enumerate(bm_names):
        event_data[f"weight_{name}"] = w_bench[:, i].astype(np.float64)

    # Nuisance variation weights (one branch per up/down variation)
    for nk, nv in nuis_weights.items():
        event_data[f"weight_{nk}"] = nv.astype(np.float64)

    # Nuisance morpher coefficients
    for i in range(a_coeffs.shape[0]):
        label = nuis_param_names[i] if i < len(nuis_param_names) else str(i)
        event_data[f"nuisance_a_{label}"] = a_coeffs[i].astype(np.float64)
        event_data[f"nuisance_b_{label}"] = b_coeffs[i].astype(np.float64)

    # Dense scale weights — one branch per (systematic, muF) point
    for syst_name, muf_dict in dense_weights.items():
        for muf_key, w in muf_dict.items():
            # muf_key is like "muf_0.5"; sanitize for ROOT branch name
            safe = muf_key.replace(".", "p")
            event_data[f"weight_{syst_name}_{safe}"] = w.astype(np.float64)

    # Sampling info
    event_data["sampling_benchmark_id"] = sampling_ids.astype(np.int32)

    # -- Metadata tree --
    meta_data = {}
    for i, name in enumerate(p_names):
        for j, bm in enumerate(bm_names):
            meta_data[f"benchmark_{bm}_{name}"] = np.array([bm_vals[j, i]])
    n_bench, n_comp = morph_mat.shape
    for i in range(n_bench):
        for j in range(n_comp):
            meta_data[f"morphing_matrix_{i}_{j}"] = np.array([morph_mat[i, j]])
    for i in range(morph_comp.shape[0]):
        for j in range(morph_comp.shape[1]):
            meta_data[f"morphing_components_{i}_{j}"] = np.array([morph_comp[i, j]], dtype=np.int32)

    with uproot.recreate(root_path) as f:
        f["Events"] = event_data
        f["Metadata"] = meta_data

    print(f"Converted to ROOT: {root_path}")
    print(f"  Events tree: {obs.shape[0]} entries, {len(event_data)} branches")
    print(f"  Metadata tree: {len(meta_data)} branches")
    print(f"  Branches:")
    for k in sorted(event_data.keys()):
        print(f"    {k}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Aggregate batch H5 files and optionally convert to ROOT")
    parser.add_argument("--input-dir", default="data/", help="Directory containing semi_parametric_samples_batch*.h5 files")
    parser.add_argument("--output", default="data/combined_samples.h5", help="Output combined HDF5 path")
    parser.add_argument("--root", action="store_true", help="Also produce a ROOT file")
    args = parser.parse_args()

    combined_path = aggregate(args.input_dir, args.output)

    if args.root:
        root_path = args.output.replace(".h5", ".root")
        convert_to_root(combined_path, root_path)
