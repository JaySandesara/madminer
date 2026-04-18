#!/usr/bin/env python
# coding: utf-8

# # Semi-Parametric Setup: End-to-End Sample Generation
# 
# This notebook generates a complete event sample with:
# - Events at the SM benchmark and all morphing basis points
# - Systematic (scale) variation weights per event (7-point grid: mu_R, mu_F, correlated)
# - Truth matrix element weights at all benchmarks (usable for morphing to arbitrary theta)
# - A standalone export file for downstream use outside MadMiner
# 
# **Batch mode:** Set `BATCH_ID` below and re-run to generate independent batches.
# Use `aggregate_and_convert.py` to combine them afterwards.

# ## 0. Environment Setup
# 
# Install the local (patched) MadMiner and set `LD_LIBRARY_PATH` so MadGraph can access LHAPDF for scale systematics. **Restart the kernel after running the pip install cell if this is your first time.**

# In[1]:


# import subprocess
# subprocess.check_call(["pip", "install", "-e", "/home/shared/madminer"])


# In[2]:


import os

os.environ["LD_LIBRARY_PATH"] = (
    "/madminer/software/MG5_aMC_v2_9_16/HEPTools/lhapdf6_py3/lib:"
    + os.environ.get("LD_LIBRARY_PATH", "")
)


# In[3]:


# ============================================================
# BATCH CONFIGURATION
# Read from env var (set by HTCondor run_semi_parametric_setup.sh) with a local default.
# ============================================================
BATCH_ID = int(os.environ.get("BATCH_ID", 1))

# ============================================================
# STAGING - large HDF5 outputs go here to avoid home quota limits
# ============================================================
STAGING_DIR = f"/staging/jsandesara/madminer/semi_parametric"
os.makedirs(STAGING_DIR, exist_ok=True)

import logging
import numpy as np
import h5py

from madminer.core import MadMiner
from madminer.lhe import LHEReader
from madminer.analysis import DataAnalyzer
from madminer.sampling import combine_and_shuffle
from particle import Particle

mg_dir = os.getenv("MG_FOLDER_PATH")

# Logging
logging.basicConfig(
    format="%(asctime)-5.5s %(name)-20.20s %(levelname)-7.7s %(message)s",
    datefmt="%H:%M",
    level=logging.INFO,
)
for key in logging.Logger.manager.loggerDict:
    if "madminer" not in key:
        logging.getLogger(key).setLevel(logging.WARNING)

print(f"Batch ID: {BATCH_ID}")


# ## 1. Parameter Space, Benchmarks, and Morphing
# 
# Define the EFT parameter space (two Wilson coefficients), manual benchmarks, and let MadMiner optimize the morphing basis.

# In[4]:


miner = MadMiner()

miner.add_parameter(
    lha_block="dim6",
    lha_id=2,
    parameter_name="CWL2",
    morphing_max_power=2,
    param_card_transform="16.52*theta",
    parameter_range=(-7.0, 7.0),
)
miner.add_parameter(
    lha_block="dim6",
    lha_id=5,
    parameter_name="CPWL2",
    morphing_max_power=2,
    param_card_transform="16.52*theta",
    parameter_range=(-7.0, 7.0),
)

# Benchmarks
miner.add_benchmark({"CWL2": 0.0, "CPWL2": 0.0}, "sm")
miner.add_benchmark({"CWL2": 5.0, "CPWL2": 0.0}, "cwl_5_cpwl_0")
miner.add_benchmark({"CWL2": -5.0, "CPWL2": 0.0}, "cwl_m5_cpwl_0")
miner.add_benchmark({"CWL2": 0.0, "CPWL2": 5.0}, "cwl_0_cpwl_5")
miner.add_benchmark({"CWL2": 0.0, "CPWL2": -5.0}, "cwl_0_cpwl_m5")
miner.add_benchmark({"CWL2": 5.0, "CPWL2": -5.0}, "cwl_5_cpwl_m5")

# Morphing (all 6 basis points already defined → fully deterministic)
miner.set_morphing(include_existing_benchmarks=True, max_overall_power=2)


# ## 2. Add Systematics and Save Setup
#
# Two nuisance parameters:
#   - scale_muf: mu_F scale variation (tree-level, no muR dependence)
#   - signal_norm: flat 10% signal normalization uncertainty
#
# PDF removed: 100 replicas blow up memory during LHE parsing and save.
# Revisit once parser is chunked or PDF is reduced to a principal-component summary.

# In[5]:


# Dense muF grid — 11 points from 0.5x to 2x nominal (only extremes get extracted by MadMiner)
miner.add_systematics(
    effect="scale",
    systematic_name="scale_muf",
    scale="muf",
    scale_variations=(0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.5, 2.0),
)

# Flat 10% signal normalization uncertainty (applied analytically, no MadGraph needed)
miner.add_systematics(
    effect="norm",
    systematic_name="signal_norm",
    norm_variation=1.10,
)

miner.save("data/setup_semi_parametric.h5")


# ## 3. Event Generation with MadGraph
# 
# Generate signal events at the SM benchmark and at each additional benchmark.
# All systematics are requested for each run so MadGraph computes the full scale variation grid.
# 
# **Note:** Change `run_card_file` to `cards/run_card_signal_large.dat` (50k events) for meaningful statistics.

# In[6]:


all_systematics = ["scale_muf", "signal_norm"]

# Only scale_muf needs MadGraph; signal_norm is applied analytically.
mg_systematics = ["scale_muf"]

additional_benchmarks = ["cwl_5_cpwl_0", "cwl_m5_cpwl_0", "cwl_0_cpwl_5", "cwl_0_cpwl_m5", "cwl_5_cpwl_m5"]
all_benchmarks = ["sm"] + additional_benchmarks

# Set to True to force regeneration even if LHE files already exist
FORCE_REGENERATE = False

def lhe_path(bm):
    return f"./mg_processes/signal_semi_{bm}_b{BATCH_ID}/Events/run_01/unweighted_events.lhe.gz"

existing = [bm for bm in all_benchmarks if os.path.exists(lhe_path(bm))]
missing = [bm for bm in all_benchmarks if bm not in existing]

if not FORCE_REGENERATE and len(missing) == 0:
    print(f"All LHE files found, skipping event generation.")
else:
    if not FORCE_REGENERATE and existing:
        print(f"LHE files found for {existing}, regenerating missing: {missing}")
        benchmarks_to_run = missing
    else:
        print(f"Running event generation for all benchmarks.")
        benchmarks_to_run = all_benchmarks

    for bm in benchmarks_to_run:
        miner.run(
            sample_benchmark=bm,
            mg_directory=mg_dir,
            mg_process_directory=f"./mg_processes/signal_semi_{bm}_b{BATCH_ID}",
            proc_card_file="cards/proc_card_signal.dat",
            param_card_template_file="cards/param_card_template.dat",
            run_card_file="cards/run_card_signal_large.dat",
            log_directory=f"logs/signal_semi_b{BATCH_ID}",
            systematics=mg_systematics,
        )


# In[8]:


# Background generation moved to semi_parametric_background.ipynb
# Run it separately — it takes much longer than signal.


# ## 4. Parse LHE Files, Define Observables, and Run Analysis
# 
# Read the generated LHE files, define detector smearing, observables, and selection cuts, then extract all event data.

# In[9]:


lhe = LHEReader("data/setup_semi_parametric.h5")

# Signal sample from SM benchmark
lhe.add_sample(
    lhe_filename=f"mg_processes/signal_semi_sm_b{BATCH_ID}/Events/run_01/unweighted_events.lhe.gz",
    sampled_from_benchmark="sm",
    is_background=False,
    k_factor=1.1,
    systematics=all_systematics,
)

# Signal samples from additional benchmarks
for bm in additional_benchmarks:
    lhe.add_sample(
        lhe_filename=f"mg_processes/signal_semi_{bm}_b{BATCH_ID}/Events/run_01/unweighted_events.lhe.gz",
        sampled_from_benchmark=bm,
        is_background=False,
        k_factor=1.1,
        systematics=all_systematics,
    )

# To include background, run semi_parametric_background.ipynb first, then uncomment:
# lhe.add_sample(
#     lhe_filename=f"mg_processes/bkg_semi_b{BATCH_ID}/Events/run_01/unweighted_events.lhe.gz",
#     sampled_from_benchmark="sm",
#     is_background=True,
#     k_factor=1.1,
#     systematics=all_systematics,
# )


# In[10]:


# Detector smearing for jet-like partons
particles = [
    *Particle.findall(lambda p: p.pdgid.is_quark),
    *Particle.findall(pdg_name="g"),
]

lhe.set_smearing(
    pdgids=[int(p.pdgid) for p in particles],
    energy_resolution_abs=0.0,
    energy_resolution_rel=0.1,
    pt_resolution_abs=None,
    pt_resolution_rel=None,
    eta_resolution_abs=0.1,
    eta_resolution_rel=0.0,
    phi_resolution_abs=0.1,
    phi_resolution_rel=0.0,
)

# Observables — full event kinematics

# Leading jet
lhe.add_observable("pt_j1",  "j[0].pt",  required=True)
lhe.add_observable("eta_j1", "j[0].eta", required=True)
lhe.add_observable("phi_j1", "j[0].phi", required=True)
lhe.add_observable("e_j1",   "j[0].e",   required=True)

# Subleading jet
lhe.add_observable("pt_j2",  "j[1].pt",  required=True)
lhe.add_observable("eta_j2", "j[1].eta", required=True)
lhe.add_observable("phi_j2", "j[1].phi", required=True)
lhe.add_observable("e_j2",   "j[1].e",   required=True)

# Dijet system
lhe.add_observable("m_jj",         "(j[0]+j[1]).m",       required=True)
lhe.add_observable("delta_eta_jj", "j[0].deltaeta(j[1])", required=True)
lhe.add_observable(
    "delta_phi_jj",
    "j[0].deltaphi(j[1]) * (-1.0 + 2.0 * float(j[0].eta > j[1].eta))",
    required=True,
)

# Leading photon
lhe.add_observable("pt_a1",  "a[0].pt",  required=True)
lhe.add_observable("eta_a1", "a[0].eta", required=True)
lhe.add_observable("phi_a1", "a[0].phi", required=True)
lhe.add_observable("e_a1",   "a[0].e",   required=True)

# Subleading photon
lhe.add_observable("pt_a2",  "a[1].pt",  required=True)
lhe.add_observable("eta_a2", "a[1].eta", required=True)
lhe.add_observable("phi_a2", "a[1].phi", required=True)
lhe.add_observable("e_a2",   "a[1].e",   required=True)

# Diphoton system
lhe.add_observable("m_aa",  "(a[0]+a[1]).m",  required=True)
lhe.add_observable("pt_aa", "(a[0]+a[1]).pt", required=True)

# Global
lhe.add_observable("met", "met.pt", required=True)

# Selection cuts
lhe.add_cut("(a[0] + a[1]).m > 122.0")
lhe.add_cut("(a[0] + a[1]).m < 128.0")
lhe.add_cut("pt_j1 > 30.0")


# In[11]:


lhe.analyse_samples()
path_to_h5 = f"lhe_data_semi_parametric_b{BATCH_ID}.h5"

path_to_saved_data = f"{STAGING_DIR}/{path_to_h5}"

lhe.save(path_to_saved_data)

# ## 5. Verify the MadMiner HDF5 File
# 
# Quick sanity check: load the file back and inspect what's inside.

# In[12]:


path_to_saved_data = f"{STAGING_DIR}/{path_to_h5}"
da = DataAnalyzer(path_to_saved_data, include_nuisance_parameters=True)

# Physical (non-nuisance) benchmark names
all_bm_names = list(da.benchmarks.keys())
bm_nuisance_flags = da.benchmark_nuisance_flags
benchmark_names_phys = [n for n, is_nuis in zip(all_bm_names, bm_nuisance_flags) if not is_nuis]

print(f"Batch {BATCH_ID}")
print(f"Parameters:          {list(da.parameters.keys())}")
print(f"Benchmarks (phys):   {benchmark_names_phys}")
print(f"Benchmarks (all):    {all_bm_names}")
print(f"Observables:         {list(da.observables.keys())}")
print(f"Systematics:         {list(da.systematics.keys())}")
print(f"Nuisance parameters: {list(da.nuisance_parameters.keys())}")
print(f"Morphing available:  {da.morpher is not None}")
print(f"Nuisance morpher:    {da.nuisance_morpher is not None}")


# In[13]:


# Load all events with benchmark weights
x_all, weights_all = da.weighted_events(theta=None)

print(f"Observations shape:  {x_all.shape}  (n_events, n_observables)")
print(f"Weights shape:       {weights_all.shape}  (n_events, n_benchmarks)")
print(f"\nFirst event observables: {x_all[0]}")
print(f"First event weights:     {weights_all[0]}")


# ## 6. Export Standalone File
# 
# Extract everything from the MadMiner HDF5 and save a self-contained file with:
# 
# | Dataset | Shape | Description |
# |---------|-------|-------------|
# | `observables` | (n_events, n_obs) | Kinematic observables per event |
# | `weights_benchmarks` | (n_events, n_benchmarks_phys) | ME weights at each physical benchmark |
# | `weights_nuisance_up` | (n_events,) | Weights at nuisance +1sigma (scale up) |
# | `weights_nuisance_down` | (n_events,) | Weights at nuisance -1sigma (scale down) |
# | `nuisance_a` | (n_nuisance, n_events) | Linear nuisance morpher coefficients |
# | `nuisance_b` | (n_nuisance, n_events) | Quadratic nuisance morpher coefficients |
# | `sampling_ids` | (n_events,) | Which benchmark each event was generated at (-1 = background) |
# | `morphing_matrix` | (n_benchmarks, n_components) | Morphing matrix for reweighting to arbitrary theta |
# | `benchmark_names` | (n_benchmarks_phys,) | Benchmark names |
# | `benchmark_values` | (n_benchmarks_phys, n_params) | Parameter values at each benchmark |
# | `parameter_names` | (n_params,) | Parameter names |
# | `observable_names` | (n_obs,) | Observable names |
# 
# With `morphing_matrix` and `benchmark_values`, you can compute `w(x|theta)` for any theta outside MadMiner.

# In[14]:


output_file = f"{STAGING_DIR}/semi_parametric_samples_batch{BATCH_ID}.h5"

# -- Collect metadata --
n_phys = len(benchmark_names_phys)
param_names = list(da.parameters.keys())
obs_names = list(da.observables.keys())

# Benchmark parameter values (n_benchmarks_phys, n_params)
benchmark_values = np.array([
    [da.benchmarks[bm].values[p] for p in param_names]
    for bm in benchmark_names_phys
])

# -- Event-level data --
x_all, w_all = da.weighted_events(theta=None)
w_phys = w_all[:, :n_phys]

# Sampling benchmark IDs
path_to_saved_data = f"{STAGING_DIR}/{path_to_h5}"
with h5py.File(path_to_saved_data, "r") as f:
    sampling_ids = f["samples/sampling_benchmarks"][()]

# -- Nuisance information --
nuisance_cols = {}
for npar_name, npar_obj in da.nuisance_parameters.items():
    if npar_obj.benchmark_pos and npar_obj.benchmark_pos in all_bm_names:
        nuisance_cols[f"{npar_name}_up"] = all_bm_names.index(npar_obj.benchmark_pos)
    if npar_obj.benchmark_neg and npar_obj.benchmark_neg in all_bm_names:
        nuisance_cols[f"{npar_name}_down"] = all_bm_names.index(npar_obj.benchmark_neg)

# Extract all nuisance weight columns
w_nuisance = {}
for col_name, col_idx in nuisance_cols.items():
    w_nuisance[col_name] = w_all[:, col_idx]

# Nuisance morpher coefficients
if da.nuisance_morpher is not None:
    nuisance_a = da.nuisance_morpher.calculate_a(w_all)
    nuisance_b = da.nuisance_morpher.calculate_b(w_all)
else:
    nuisance_a = np.zeros((0, x_all.shape[0]))
    nuisance_b = np.zeros((0, x_all.shape[0]))

# -- Morphing matrix --
morphing_matrix = da.morpher.morphing_matrix if da.morpher is not None else None
morphing_components = da.morpher.components if da.morpher is not None else None

print(f"Batch {BATCH_ID}")
print(f"Events:              {x_all.shape[0]}")
print(f"Observables:         {x_all.shape[1]} {obs_names}")
print(f"Physical benchmarks: {n_phys} {benchmark_names_phys}")
print(f"Nuisance columns:    {list(nuisance_cols.keys())}")
print(f"Nuisance a shape:    {nuisance_a.shape}")
print(f"Morphing matrix:     {morphing_matrix.shape if morphing_matrix is not None else 'N/A'}")


# In[15]:


# -- Write standalone HDF5 --
with h5py.File(output_file, "w") as f:
    # Event data
    f.create_dataset("observables", data=x_all)
    f.create_dataset("weights_benchmarks", data=w_phys)
    f.create_dataset("sampling_ids", data=sampling_ids)

    # Nuisance weights (one dataset per variation)
    nuis_grp = f.create_group("nuisance_weights")
    for col_name, col_data in w_nuisance.items():
        nuis_grp.create_dataset(col_name, data=col_data)

    # Nuisance morpher coefficients
    f.create_dataset("nuisance_a", data=nuisance_a)
    f.create_dataset("nuisance_b", data=nuisance_b)

    # Morphing
    if morphing_matrix is not None:
        f.create_dataset("morphing_matrix", data=morphing_matrix)
    if morphing_components is not None:
        f.create_dataset("morphing_components", data=morphing_components)

    # Metadata
    f.create_dataset("benchmark_names", data=np.array(benchmark_names_phys, dtype="S256"))
    f.create_dataset("benchmark_values", data=benchmark_values)
    f.create_dataset("parameter_names", data=np.array(param_names, dtype="S256"))
    f.create_dataset("observable_names", data=np.array(obs_names, dtype="S256"))
    f.create_dataset("nuisance_parameter_names", data=np.array(list(da.nuisance_parameters.keys()), dtype="S256"))
    f.create_dataset("systematics_names", data=np.array(list(da.systematics.keys()), dtype="S256"))
    f.attrs["batch_id"] = BATCH_ID

    # Dense scale weights — raw per-event muF variation weights.
    # Apply the same sampling_factor that DataAnalyzer applies to benchmark weights
    # (see DataAnalyzer._calculate_sampling_factors + load_events in hdf5.py):
    # each event's weight is multiplied by n_events_generated_at_its_benchmark / n_total.
    with h5py.File(path_to_saved_data, "r") as src:
        if "scale_weights_dense" in src:
            events_per_bm = np.asarray(da.n_events_generated_per_benchmark, dtype=np.float64)
            sampling_factors = events_per_bm / events_per_bm.sum()
            sampling_factors = np.hstack((sampling_factors, 1.0))  # background slot
            per_event_factor = sampling_factors[sampling_ids]

            dgrp = f.create_group("scale_weights_dense")
            for syst_name in src["scale_weights_dense"]:
                ssrc = src[f"scale_weights_dense/{syst_name}"]
                sdst = dgrp.create_group(syst_name)
                for key in ssrc:
                    sdst.create_dataset(key, data=ssrc[key][()] * per_event_factor)
            print(f"Included dense scale weights: {list(src['scale_weights_dense'].keys())}")

print(f"Saved batch {BATCH_ID} to {output_file}")


# ## 7. Verify: Reconstruct Weights at Arbitrary Theta
# 
# Demonstrate that the exported file is self-contained: load it back without MadMiner and compute `w(x|theta)` at an arbitrary parameter point using the morphing matrix.

# In[16]:


# Load the standalone file and verify
with h5py.File(output_file, "r") as f:
    print(f"Batch {f.attrs['batch_id']}")
    print(f"Events: {f['observables'].shape[0]}")
    print(f"Observables: {f['observables'].shape[1]}")
    print(f"Benchmark weights: {f['weights_benchmarks'].shape}")
    print(f"Nuisance weight datasets: {list(f['nuisance_weights'].keys())}")
    print(f"Nuisance a/b shape: {f['nuisance_a'].shape}")
    print(f"Morphing matrix: {f['morphing_matrix'].shape}")
    print()

    # Quick validation: check all weights are nonzero
    w = f['weights_benchmarks'][()]
    bm_names = [s.decode() for s in f['benchmark_names'][()]]
    for i, name in enumerate(bm_names):
        col = w[:, i]
        print(f"  {name:30s}  mean={col.mean():.4e}  nonzero={np.count_nonzero(col)}/{len(col)}")
    print()
    for nk in f['nuisance_weights'].keys():
        nw = f['nuisance_weights'][nk][()]
        print(f"  {nk:40s}  mean={nw.mean():.4e}  nonzero={np.count_nonzero(nw)}/{len(nw)}")


# In[ ]:





# In[ ]:





# In[ ]:




