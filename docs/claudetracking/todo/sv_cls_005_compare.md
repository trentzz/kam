# SV-CLS-005: Compare Before and After Classification Fix

**Epic**: SV-CLASSIFY
**Priority**: medium
**Depends on**: SV-CLS-004, SV-BASE-007
**Status**: todo

## Goal

Score the re-run NovelInsertion results and compare sensitivity before and
after the classification fix.

## Success Criteria

- [ ] Comparison table showing classification counts and sensitivity before/after the fix.
- [ ] Results documented in the investigation document from SV-CLS-001.
- [ ] `/update` has been run after changes.

## Steps

1. Score the re-run results using score_sv_new_suite.py.
2. Compare against the original (pre-fix) scored results.
3. Produce a summary table.

## Notes

Sensitivity numbers may or may not change, since the underlying detection
may have been working but labelled incorrectly. The key metric is correct
SVTYPE assignment.
