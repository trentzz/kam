# SV-SWEEP: Parameter and Size Sweep Experiments

**Status**: todo
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
| SV-SWP-001 | todo/sv_swp_001_kmer_script.md | todo |
| SV-SWP-002 | todo/sv_swp_002_kmer_run.md | todo |
| SV-SWP-003 | todo/sv_swp_003_threshold_script.md | todo |
| SV-SWP-004 | todo/sv_swp_004_threshold_run.md | todo |
| SV-SWP-005 | todo/sv_swp_005_size_configs.md | todo |
| SV-SWP-006 | todo/sv_swp_006_size_run.md | todo |
| SV-SWP-007 | todo/sv_swp_007_aggregate.md | todo |

## Scope

- `docs/benchmarking/sv_new/scripts/` — sweep runner and scoring scripts
- `docs/benchmarking/sv_new/configs/` — size-sweep configs
- `docs/benchmarking/sv_new/results/sweep_*/` — sweep results
- `docs/benchmarking/sv_new/figures/` — heatmaps and line plots
