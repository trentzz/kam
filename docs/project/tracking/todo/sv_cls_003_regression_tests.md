# SV-CLS-003: Add Regression Tests

**Epic**: SV-CLASSIFY
**Priority**: high
**Depends on**: SV-CLS-002
**Status**: todo

## Goal

Add at least 3 regression tests using real benchmark sequences to prevent
future misclassification regressions.

## Success Criteria

- [ ] Test in kam-call covering the exact benchmark NovelInsertion case.
- [ ] Test covering a true TandemDuplication case that must remain classified as DUP.
- [ ] Test covering a borderline case near the classification boundary.
- [ ] All tests pass with `cargo test`.
- [ ] `/update` has been run after changes.

## Steps

1. Extract representative sequences from benchmark results for each case.
2. Write test functions in the appropriate test module in kam-call.
3. Verify all tests pass.

## Notes

Use sequences from actual varforge-generated data to ensure tests are
realistic rather than synthetic edge cases.
