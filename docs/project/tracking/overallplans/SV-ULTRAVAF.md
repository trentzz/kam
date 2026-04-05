# SV-ULTRAVAF: Ultra-Low VAF Detection Limits

**Status**: todo
**Priority**: medium-low
**Branch**: epic/SV-ULTRAVAF

## Goal

Extend the VAF range below 0.05% with ultra-low VAF configs (0.01%, 0.02%,
0.03%, 0.04%) to establish the detection floor for each SV type.

## Motivation

The current benchmark suite starts at 0.05% VAF. For clinical applications
(MRD monitoring), detection at 0.01-0.04% VAF is relevant. Knowing where
sensitivity drops to zero guides protocol design.

## Child tasks

| ID | File | Status |
|----|------|--------|
| SV-ULV-001 | todo/sv_ulv_001_configs.md | todo |
| SV-ULV-002 | todo/sv_ulv_002_sims.md | todo |
| SV-ULV-003 | todo/sv_ulv_003_kam_score.md | todo |
| SV-ULV-004 | todo/sv_ulv_004_curves.md | todo |

## Scope

- `docs/benchmarking/sv_new/configs/` — ultra-low VAF configs
- `docs/benchmarking/sv_new/results/` — ultra-low VAF results
- `docs/benchmarking/sv_new/figures/` — detection limit curves
