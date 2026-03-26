# SV-BASE-007: Score All Baseline Results

**Epic**: SV-BASELINE
**Priority**: critical
**Depends on**: SV-BASE-004, SV-BASE-005, SV-BASE-006
**Status**: todo

## Goal

Run score_sv_new_suite.py on all InvDel, NovelInsertion, and Fusion results to
produce per-dataset and summary scoring tables.

## Success Criteria

- [ ] sv_new_per_dataset.tsv exists in `docs/benchmarking/sv_new/summary/` with rows for all 150 datasets.
- [ ] sv_new_summary.tsv exists with sensitivity at key VAFs per type and mode.
- [ ] `/update` has been run after changes.

## Steps

1. Run score_sv_new_suite.py against all result directories.
2. Inspect output for completeness and sanity.
3. Record any anomalies (e.g. zero sensitivity at high VAFs).

## Notes

The scoring script compares kam VCF calls against truth VCFs using position
and type matching. Check that the matching logic handles all three SV types.
