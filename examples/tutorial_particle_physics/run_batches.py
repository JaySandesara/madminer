"""
Run event generation and export for multiple batches.

Usage:
    python run_batches.py --start 1 --end 5
    python run_batches.py --start 1 --end 5 --run-card cards/run_card_signal_large.dat
    python run_batches.py --start 3 --end 3  # single batch

Runs the full pipeline per batch:
  1. MadGraph event generation (signal at all benchmarks + background)
  2. LHE parsing with observables, cuts, smearing
  3. Export to standalone HDF5

After all batches complete, run:
    python aggregate_and_convert.py --input-dir data/ --output data/combined_samples.h5 --root
"""

import argparse
import logging
import os
import sys

import numpy as np
import h5py

from madminer.core import MadMiner
from madminer.lhe import LHEReader
from madminer.analysis import DataAnalyzer
from particle import Particle


def setup_logging():
    logging.basicConfig(
        format="%(asctime)-5.5s %(name)-20.20s %(levelname)-7.7s %(message)s",
        datefmt="%H:%M",
        level=logging.INFO,
    )
    for key in logging.Logger.manager.loggerDict:
        if "madminer" not in key:
            logging.getLogger(key).setLevel(logging.WARNING)


def create_setup():
    """Create parameter space, benchmarks, morphing, and systematics. Only needs to run once."""
    setup_file = "data/setup_semi_parametric.h5"

    miner = MadMiner()

    miner.add_parameter(
        lha_block="dim6", lha_id=2, parameter_name="CWL2",
        morphing_max_power=2, param_card_transform="16.52*theta",
        parameter_range=(-20.0, 20.0),
    )
    miner.add_parameter(
        lha_block="dim6", lha_id=5, parameter_name="CPWL2",
        morphing_max_power=2, param_card_transform="16.52*theta",
        parameter_range=(-20.0, 20.0),
    )

    miner.add_benchmark({"CWL2": 0.0, "CPWL2": 0.0}, "sm")
    miner.add_benchmark({"CWL2": 15.2, "CPWL2": 0.1}, "w")
    miner.add_benchmark({"CWL2": -15.4, "CPWL2": 0.2}, "neg_w")
    miner.add_benchmark({"CWL2": 0.3, "CPWL2": 15.1}, "ww")
    miner.add_benchmark({"CWL2": 0.4, "CPWL2": -15.3}, "neg_ww")

    miner.set_morphing(include_existing_benchmarks=True, max_overall_power=2)

    miner.add_systematics(effect="scale", systematic_name="scale_mur", scale="mur")
    miner.add_systematics(effect="scale", systematic_name="scale_muf", scale="muf")
    miner.add_systematics(effect="scale", systematic_name="scale_corr", scale="mu")

    miner.save(setup_file)
    return miner, setup_file


def run_generation(miner, batch_id, signal_run_card, bkg_run_card, skip_background):
    """Run MadGraph event generation for one batch."""
    mg_dir = os.getenv("MG_FOLDER_PATH")
    if not mg_dir:
        print("ERROR: MG_FOLDER_PATH not set")
        sys.exit(1)

    # Only scale_corr needed for MadGraph (it computes full 3x3 grid)
    mg_systematics = ["scale_corr"]

    # Signal at SM
    print(f"\n[Batch {batch_id}] Generating signal at SM...")
    miner.run(
        sample_benchmark="sm",
        mg_directory=mg_dir,
        mg_process_directory=f"./mg_processes/signal_semi_sm_b{batch_id}",
        proc_card_file="cards/proc_card_signal.dat",
        param_card_template_file="cards/param_card_template.dat",
        run_card_file=signal_run_card,
        log_directory=f"logs/signal_semi_b{batch_id}",
        systematics=mg_systematics,
    )

    # Signal at additional benchmarks
    additional_benchmarks = ["w", "ww", "neg_w", "neg_ww"]
    for bm in additional_benchmarks:
        print(f"\n[Batch {batch_id}] Generating signal at {bm}...")
        miner.run(
            sample_benchmark=bm,
            mg_directory=mg_dir,
            mg_process_directory=f"./mg_processes/signal_semi_{bm}_b{batch_id}",
            proc_card_file="cards/proc_card_signal.dat",
            param_card_template_file="cards/param_card_template.dat",
            run_card_file=signal_run_card,
            log_directory=f"logs/signal_semi_b{batch_id}",
            systematics=mg_systematics,
        )

    # Background
    if not skip_background:
        print(f"\n[Batch {batch_id}] Generating background...")
        miner.run(
            sample_benchmark="sm",
            mg_directory=mg_dir,
            mg_process_directory=f"./mg_processes/bkg_semi_b{batch_id}",
            proc_card_file="cards/proc_card_background.dat",
            param_card_template_file="cards/param_card_template.dat",
            run_card_file=bkg_run_card,
            log_directory=f"logs/bkg_semi_b{batch_id}",
            systematics=mg_systematics,
        )

    return additional_benchmarks


def run_analysis(batch_id, additional_benchmarks, setup_file, skip_background):
    """Parse LHE files and extract observables + weights."""
    all_systematics = ["scale_mur", "scale_muf", "scale_corr"]

    lhe = LHEReader(setup_file)

    # Signal at SM
    lhe.add_sample(
        lhe_filename=f"mg_processes/signal_semi_sm_b{batch_id}/Events/run_01/unweighted_events.lhe.gz",
        sampled_from_benchmark="sm",
        is_background=False,
        k_factor=1.1,
        systematics=all_systematics,
    )

    # Signal at additional benchmarks
    for bm in additional_benchmarks:
        lhe.add_sample(
            lhe_filename=f"mg_processes/signal_semi_{bm}_b{batch_id}/Events/run_01/unweighted_events.lhe.gz",
            sampled_from_benchmark=bm,
            is_background=False,
            k_factor=1.1,
            systematics=all_systematics,
        )

    # Background
    if not skip_background:
        lhe.add_sample(
            lhe_filename=f"mg_processes/bkg_semi_b{batch_id}/Events/run_01/unweighted_events.lhe.gz",
            sampled_from_benchmark="sm",
            is_background=True,
            k_factor=1.1,
            systematics=all_systematics,
        )

    # Smearing
    particles = [
        *Particle.findall(lambda p: p.pdgid.is_quark),
        *Particle.findall(pdg_name="g"),
    ]
    lhe.set_smearing(
        pdgids=[int(p.pdgid) for p in particles],
        energy_resolution_abs=0.0, energy_resolution_rel=0.1,
        pt_resolution_abs=None, pt_resolution_rel=None,
        eta_resolution_abs=0.1, eta_resolution_rel=0.0,
        phi_resolution_abs=0.1, phi_resolution_rel=0.0,
    )

    # Observables
    lhe.add_observable("pt_j1", "j[0].pt", required=False, default=0.0)
    lhe.add_observable(
        "delta_phi_jj",
        "j[0].deltaphi(j[1]) * (-1.0 + 2.0 * float(j[0].eta > j[1].eta))",
        required=True,
    )
    lhe.add_observable("met", "met.pt", required=True)

    # Cuts
    lhe.add_cut("(a[0] + a[1]).m > 122.0")
    lhe.add_cut("(a[0] + a[1]).m < 128.0")
    lhe.add_cut("pt_j1 > 30.0")

    madminer_file = f"data/lhe_data_semi_parametric_b{batch_id}.h5"
    lhe.analyse_samples()
    lhe.save(madminer_file)

    return madminer_file


def export_standalone(batch_id, madminer_file):
    """Export standalone HDF5 with all weights and metadata."""
    da = DataAnalyzer(madminer_file, include_nuisance_parameters=True)

    all_bm_names = list(da.benchmarks.keys())
    bm_nuisance_flags = da.benchmark_nuisance_flags
    benchmark_names_phys = [n for n, is_nuis in zip(all_bm_names, bm_nuisance_flags) if not is_nuis]
    n_phys = len(benchmark_names_phys)
    param_names = list(da.parameters.keys())
    obs_names = list(da.observables.keys())

    benchmark_values = np.array([
        [da.benchmarks[bm].values[p] for p in param_names]
        for bm in benchmark_names_phys
    ])

    x_all, w_all = da.weighted_events(theta=None)
    w_phys = w_all[:, :n_phys]

    with h5py.File(madminer_file, "r") as f:
        sampling_ids = f["samples/sampling_benchmarks"][()]

    # Nuisance weights
    nuisance_cols = {}
    for npar_name, npar_obj in da.nuisance_parameters.items():
        if npar_obj.benchmark_pos and npar_obj.benchmark_pos in all_bm_names:
            nuisance_cols[f"{npar_name}_up"] = all_bm_names.index(npar_obj.benchmark_pos)
        if npar_obj.benchmark_neg and npar_obj.benchmark_neg in all_bm_names:
            nuisance_cols[f"{npar_name}_down"] = all_bm_names.index(npar_obj.benchmark_neg)

    w_nuisance = {col_name: w_all[:, col_idx] for col_name, col_idx in nuisance_cols.items()}

    # Morpher coefficients
    if da.nuisance_morpher is not None:
        nuisance_a = da.nuisance_morpher.calculate_a(w_all)
        nuisance_b = da.nuisance_morpher.calculate_b(w_all)
    else:
        nuisance_a = np.zeros((0, x_all.shape[0]))
        nuisance_b = np.zeros((0, x_all.shape[0]))

    morphing_matrix = da.morpher.morphing_matrix if da.morpher is not None else None
    morphing_components = da.morpher.components if da.morpher is not None else None

    # Write
    output_file = f"data/semi_parametric_samples_batch{batch_id}.h5"
    with h5py.File(output_file, "w") as f:
        f.create_dataset("observables", data=x_all)
        f.create_dataset("weights_benchmarks", data=w_phys)
        f.create_dataset("sampling_ids", data=sampling_ids)

        nuis_grp = f.create_group("nuisance_weights")
        for col_name, col_data in w_nuisance.items():
            nuis_grp.create_dataset(col_name, data=col_data)

        f.create_dataset("nuisance_a", data=nuisance_a)
        f.create_dataset("nuisance_b", data=nuisance_b)

        if morphing_matrix is not None:
            f.create_dataset("morphing_matrix", data=morphing_matrix)
        if morphing_components is not None:
            f.create_dataset("morphing_components", data=morphing_components)

        f.create_dataset("benchmark_names", data=np.array(benchmark_names_phys, dtype="S256"))
        f.create_dataset("benchmark_values", data=benchmark_values)
        f.create_dataset("parameter_names", data=np.array(param_names, dtype="S256"))
        f.create_dataset("observable_names", data=np.array(obs_names, dtype="S256"))
        f.create_dataset("nuisance_parameter_names", data=np.array(list(da.nuisance_parameters.keys()), dtype="S256"))
        f.create_dataset("systematics_names", data=np.array(list(da.systematics.keys()), dtype="S256"))
        f.attrs["batch_id"] = batch_id

    print(f"[Batch {batch_id}] Exported {x_all.shape[0]} events to {output_file}")
    return output_file


def run_batch(batch_id, signal_run_card, bkg_run_card, skip_background):
    """Run the full pipeline for a single batch."""
    print(f"\n{'='*60}")
    print(f"  BATCH {batch_id}")
    print(f"{'='*60}")

    miner, setup_file = create_setup()
    additional_benchmarks = run_generation(miner, batch_id, signal_run_card, bkg_run_card, skip_background)
    madminer_file = run_analysis(batch_id, additional_benchmarks, setup_file, skip_background)
    output_file = export_standalone(batch_id, madminer_file)

    return output_file


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run MadMiner batch generation")
    parser.add_argument("--start", type=int, default=1, help="First batch ID")
    parser.add_argument("--end", type=int, default=1, help="Last batch ID (inclusive)")
    parser.add_argument("--signal-run-card", default="cards/run_card_signal_large.dat")
    parser.add_argument("--bkg-run-card", default="cards/run_card_background.dat")
    parser.add_argument("--skip-background", action="store_true", help="Skip background generation")
    args = parser.parse_args()

    # Environment
    os.environ["LD_LIBRARY_PATH"] = (
        "/madminer/software/MG5_aMC_v2_9_16/HEPTools/lhapdf6_py3/lib:"
        + os.environ.get("LD_LIBRARY_PATH", "")
    )

    setup_logging()

    output_files = []
    for batch_id in range(args.start, args.end + 1):
        out = run_batch(batch_id, args.signal_run_card, args.bkg_run_card, args.skip_background)
        output_files.append(out)

    print(f"\n{'='*60}")
    print(f"  ALL DONE — {len(output_files)} batches")
    print(f"{'='*60}")
    for f in output_files:
        print(f"  {f}")
    print(f"\nNext: python aggregate_and_convert.py --input-dir data/ --output data/combined_samples.h5 --root")
