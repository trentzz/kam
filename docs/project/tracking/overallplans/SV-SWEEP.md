# SV-SWEEP: Parameter and Size Sweep Experiments

**Status**: in-progress
**Priority**: medium
**Branch**: epic/SV-SWEEP

## Goal

Systematically explore the parameter space (k-mer size, confidence threshold,
min alt molecules) and SV size space (20-500 bp) to map detection boundaries
and find optimal configurations per SV type.

## Motivation

The current defaults (k=31, sv_min_confidence=0.95, sv_min_alt_molecules=1)
were chosen by reasoning, not by empirical optimisation. Sweeping these
parameters on synthetic data quantifies the tradeoffs and may reveal better
operating points.

## Child tasks

| ID | File | Status |
|----|------|--------|
| SV-SWP-001 | done/sv_swp_001_kmer_script.md | done |
| SV-SWP-002 | todo/sv_swp_002_kmer_run.md | todo |
| SV-SWP-003 | done/sv_swp_003_threshold_script.md | done |
| SV-SWP-004 | todo/sv_swp_004_threshold_run.md | todo |
| SV-SWP-005 | done/sv_swp_005_size_configs.md | done |
| SV-SWP-006 | todo/sv_swp_006_size_run.md | todo |
| SV-SWP-007 | done/sv_swp_007_aggregate.md | done |

## Scripts

| Script | Purpose |
|--------|---------|
| `scripts/benchmarking/sv_kmer_sweep.sh` | K-mer size sweep (k=21,25,27,31,35,41) |
| `scripts/benchmarking/sv_threshold_sweep.sh` | Confidence and min-alt-molecules sweep |
| `scripts/benchmarking/sv_detection_limits.sh` | Ultra-low VAF detection limit runner |
| `scripts/benchmarking/aggregate_sweep.py` | Aggregate results, produce heatmaps and recommendations |
| `scripts/benchmarking/generate_size_sweep_configs.py` | Generate size sweep varforge configs |
| `scripts/benchmarking/generate_ultra_low_vaf_configs.py` | Generate ultra-low VAF varforge configs |

## Configs

| Directory | Contents |
|-----------|----------|
| `scripts/benchmarking/sv_size_sweep_configs/` | 28 configs: DEL/DUP at 5 sizes, INV at 4 sizes (>= 50 bp), 2 replicates each |
| `docs/benchmarking/03-sv-extended/configs/ultra_low_vaf/` | 48 configs: 6 SV types x 4 VAF levels (0.01-0.04%) x 2 replicates |

## Scope

- `scripts/benchmarking/` — sweep runner, scoring, and aggregation scripts
- `scripts/benchmarking/sv_size_sweep_configs/` — size-sweep varforge configs
- `docs/benchmarking/03-sv-extended/configs/ultra_low_vaf/` — ultra-low VAF configs
- `docs/benchmarking/03-sv-extended/` — extended SV benchmark data and results
