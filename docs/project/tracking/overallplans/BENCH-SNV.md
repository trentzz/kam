# BENCH-SNV: SNV and Indel Benchmarking

**Status**: done
**Priority**: high

## Goal

Build a complete benchmarking suite for SNV and indel detection in kam. This
covers:

1. Separating SNV and indel metrics in the scoring pipeline.
2. Synthetic datasets generated with varforge for unit-level benchmarking
   (small, self-contained, no real patient data).
3. Configs for the synthetic datasets committed to the repo; the generated
   FASTQs are gitignored.

## Motivation

The current titration benchmark reports combined sensitivity. SNV and indel
sensitivity differ substantially (SNVs ~70-80%, indels ~31-39% at 2% VAF).
Reporting them separately gives a clearer picture of where kam is strong and
where it needs work. Indels are a known weakness due to the anchor coverage
problem.

Synthetic datasets let us run fast regression tests without needing the full
Twist titration FASTQs, which are large and not in the repo.

## Child tasks

| ID | File | Status |
|----|------|--------|
| BENCH-SNV-001 | done/bench_snv_001_separate_scoring.md | done |
| BENCH-SNV-002 | done/bench_snv_002_synthetic_datasets.md | done |
| BENCH-SNV-003 | done/bench_snv_003_snv_vaf_sweep.md | done |
| BENCH-SNV-004 | done/bench_snv_004_indel_vaf_sweep.md | done |
| BENCH-SNV-005 | done/bench_snv_005_combined_datasets.md | done |

## Scope

- `docs/benchmarking/snvindel/scripts/score_variants.py` — SNV/indel split
- `docs/benchmarking/snvindel/scripts/configs/` — varforge configs for synthetic data
- `docs/benchmarking/snvindel/data/` — small synthetic reference and truth VCFs

## Out of scope

- Changes to the kam Rust code
- Re-running the full titration benchmark
