# SV-BASELINE: Complete Baseline SV Benchmark Suite

**Status**: todo
**Priority**: critical
**Branch**: epic/SV-BASELINE

## Goal

Finish all remaining varforge simulations and kam runs for InvDel, NovelInsertion,
and Fusion. Score everything and produce the first complete baseline numbers with
sensitivity curves.

## Motivation

SV-EXPAND added classification logic and BENCH-SV-NEW generated configs, but the
actual benchmark runs have not been executed. Without baseline numbers, no further
experiments (parameter sweeps, code fixes) have a point of comparison.

## Child tasks

| ID | File | Status |
|----|------|--------|
| SV-BASE-001 | todo/sv_base_001_invdel_sims.md | todo |
| SV-BASE-002 | todo/sv_base_002_novins_sims.md | todo |
| SV-BASE-003 | todo/sv_base_003_fusion_sims.md | todo |
| SV-BASE-004 | todo/sv_base_004_invdel_kam.md | todo |
| SV-BASE-005 | todo/sv_base_005_novins_kam.md | todo |
| SV-BASE-006 | todo/sv_base_006_fusion_kam.md | todo |
| SV-BASE-007 | todo/sv_base_007_score.md | todo |
| SV-BASE-008 | todo/sv_base_008_figures.md | todo |

## Scope

- `docs/benchmarking/sv_new/results/` — varforge sim output and kam results
- `docs/benchmarking/sv_new/summary/` — scoring output
- `docs/benchmarking/sv_new/figures/` — sensitivity curves
