"""
Convert semi_parametric_samples.h5 to ROOT format.

Usage:
    python convert_to_root.py [--input data/semi_parametric_samples.h5] [--output data/semi_parametric_samples.root]

Requires: uproot, h5py, numpy
    pip install uproot h5py numpy
"""

import argparse
import h5py
import numpy as np

try:
    import uproot
except ImportError:
    raise ImportError("uproot is required: pip install uproot")


def convert(input_path, output_path):
    with h5py.File(input_path, "r") as f:
        obs = f["observables"][()]
        w_bench = f["weights_benchmarks"][()]
        w_up = f["weights_nuisance_up"][()]
        w_down = f["weights_nuisance_down"][()]
        a_coeffs = f["nuisance_a"][()]
        b_coeffs = f["nuisance_b"][()]
        sampling_ids = f["sampling_ids"][()]
        morph_mat = f["morphing_matrix"][()]
        morph_comp = f["morphing_components"][()]
        bm_names = [s.decode() for s in f["benchmark_names"][()]]
        bm_vals = f["benchmark_values"][()]
        p_names = [s.decode() for s in f["parameter_names"][()]]
        o_names = [s.decode() for s in f["observable_names"][()]]

    # Build event tree
    event_data = {}

    # Observables
    for i, name in enumerate(o_names):
        event_data[name] = obs[:, i].astype(np.float64)

    # ME weights at each physical benchmark
    for i, name in enumerate(bm_names):
        event_data[f"weight_{name}"] = w_bench[:, i].astype(np.float64)

    # Nuisance weights
    event_data["weight_scale_up"] = w_up.astype(np.float64)
    event_data["weight_scale_down"] = w_down.astype(np.float64)

    # Nuisance morpher coefficients
    for i in range(a_coeffs.shape[0]):
        event_data[f"nuisance_a_{i}"] = a_coeffs[i].astype(np.float64)
        event_data[f"nuisance_b_{i}"] = b_coeffs[i].astype(np.float64)

    # Sampling info
    event_data["sampling_benchmark_id"] = sampling_ids.astype(np.int32)

    # Build metadata tree
    meta_data = {}
    for i, name in enumerate(p_names):
        for j, bm in enumerate(bm_names):
            meta_data[f"benchmark_{bm}_{name}"] = np.array([bm_vals[j, i]])

    # Morphing matrix (flattened into branches)
    n_bench, n_comp = morph_mat.shape
    for i in range(n_bench):
        for j in range(n_comp):
            meta_data[f"morphing_matrix_{i}_{j}"] = np.array([morph_mat[i, j]])

    for i in range(morph_comp.shape[0]):
        for j in range(morph_comp.shape[1]):
            meta_data[f"morphing_components_{i}_{j}"] = np.array([morph_comp[i, j]], dtype=np.int32)

    with uproot.recreate(output_path) as f:
        f["Events"] = event_data
        f["Metadata"] = meta_data

    n_events = obs.shape[0]
    print(f"Converted {n_events} events to {output_path}")
    print(f"  Events tree: {len(event_data)} branches")
    print(f"    Observables:  {o_names}")
    print(f"    ME weights:   weight_{{{'|'.join(bm_names)}}}")
    print(f"    Scale:        weight_scale_up, weight_scale_down")
    print(f"    Nuisance:     nuisance_a_*, nuisance_b_*")
    print(f"    Sampling:     sampling_benchmark_id")
    print(f"  Metadata tree: {len(meta_data)} branches (benchmark values + morphing matrix)")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert semi_parametric_samples.h5 to ROOT")
    parser.add_argument("--input", default="data/semi_parametric_samples.h5")
    parser.add_argument("--output", default="data/semi_parametric_samples.root")
    args = parser.parse_args()
    convert(args.input, args.output)