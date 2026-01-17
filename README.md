# graphSSR: SSR detection and parameter optimization toolkit

This directory contains GraphSSR (C++/Python) and a set of parameter-optimization
scripts for SSR callers (TRF, Ultra, PyTRF, and custom finders).

## Contents

- `graphSSR` (binary): C++ SSR finder
- `graphSSR.cpp`: C++ source for GraphSSR
- `graphSSR.py`: Python wrapper/runner for GraphSSR
- `optimize_graph.py`: parameter optimization for GraphSSR
- `optimize_trf.py`: parameter optimization for TRF
- `optimize_ultra.py`: parameter optimization for Ultra
- `optimize_pytrf.py`: parameter optimization for PyTRF
- `optimize_finders.py`: multi-finder optimization/benchmarking

## Requirements

- Linux
- Python 3.8+ (recommended)
- TRF (Tandem Repeats Finder) for `optimize_trf.py`
- PyTRF for `optimize_pytrf.py`
- A C++17 compiler (if you want to rebuild `graphSSR`)
- Common Python packages (depending on scripts): `numpy`, `pandas`

## Build GraphSSR (optional)

```bash
g++ -O3 -std=c++17 -o graphSSR graphSSR.cpp
```

## Example usage

### TRF parameter optimization

```bash
python3 optimize_trf.py \
  --trf_path /path/to/trf409.linux64 \
  --alignment_sets "2,5,5;2,7,7;2,9,9" \
  --pm_range "75,80" \
  --pi_range "10,15,20" \
  --score_range "30,40,50" \
  --num_replicates 1000 \
  --max_workers 4 \
  --seed 1 \
  --match_mode interval \
  --reuse_sim \
  --output_csv trf_grid_search_results.csv \
  --calibrate_threshold
```

### GraphSSR parameter optimization

```bash
python3 optimize_graph.py \
  --graph_path /path/to/graphSSR \
  --graph_k_range "5,7,9,11,13,15" \
  --graph_min_len_range "10,20,30" \
  --graph_min_cov_range "1.0,2.0,3.0,4.0,5.0" \
  --graph_max_merge_dist_range "15,20,25" \
  --num_replicates 1000 \
  --max_workers 4 \
  --seed 1 \
  --reuse_sim \
  --match_mode interval \
  --calibrate_threshold \
  --output_csv graph_grid_search_results.csv
```

### Ultra parameter optimization

```bash
python3 optimize_ultra.py \
  --match_range "0.7,0.8,0.9" \
  --decay_range "0.8,0.85,0.9" \
  --ri_range "0.01,0.02,0.03" \
  --rd_range "0.01,0.02,0.03" \
  --num_replicates 1000 \
  --max_workers 4 \
  --seed 1 \
  --match_mode interval \
  --reuse_sim \
  --calibrate_threshold \
  --output_csv ultra_grid_search_results.csv
```

### PyTRF parameter optimization

```bash
python3 optimize_pytrf.py \
  --min_identity_range "0.7,0.8,0.9" \
  --max_error_range "1,2,3,4" \
  --min_seed_repeat_range "2,3,4,5" \
  --num_replicates 1000 \
  --max_workers 4 \
  --seed 1 \
  --match_mode interval \
  --reuse_sim \
  --calibrate_threshold \
  --output_csv pytrf_grid_search_results.csv
```

## Outputs

Most optimization scripts produce:
- `*_grid_search_results.csv`
- optional threshold selection metadata/curves

## Notes

- Paths to external tools must be set explicitly (e.g., `--trf_path`).
- Match modes (e.g., `interval` vs. `interval+period`) should be reported in methods.

## License

CC BY 4.0

## Citation

In preparation.
