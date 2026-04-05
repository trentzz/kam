# BENCH-VF-005: Score and Aggregate Results

**Epic**: BENCH-VARFORGE (docs/claudetracking/overallplans/BENCH-VARFORGE.md)
**Priority**: high
**Depends on**: BENCH-VF-004
**Status**: todo

## Goal

Score all 50-dataset runs per variant class. Produce per-dataset sensitivity
tables and aggregate summary tables for discovery and tumour-informed modes.

## Steps

1. Write `docs/benchmarking/snvindel/scripts/score_snv_indel_suite.py` to:
   - Score each dataset (TP/FP/FN, sensitivity, precision) for both modes.
   - Write a per-dataset TSV with columns: type, vaf, replicate, mode,
     tp, fp, fn, sensitivity, precision, f1.
   - Write a summary TSV with sensitivity at key VAF levels (0.1%, 0.25%,
     0.5%, 1%, 2%, 5%) for each mode.
2. Write an equivalent `score_sv_suite.py` for SV datasets.
3. Run both scripts and commit the output TSVs.

## Notes

Use the existing `score_variants.py` as the per-dataset scoring function.
Matching is by (REF, ALT) sequence only.
Both TSVs must be committed to the repo (they are small text files, not
gitignored). Place in `docs/benchmarking/snvindel/summary/` and
`docs/benchmarking/sv/summary/`.
