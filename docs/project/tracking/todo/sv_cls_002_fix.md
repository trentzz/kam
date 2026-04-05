# SV-CLS-002: Fix Classification Logic

**Epic**: SV-CLASSIFY
**Priority**: high
**Depends on**: SV-CLS-001
**Status**: todo

## Goal

Fix `is_tandem_duplication()` or surrounding logic so that novel insertions
are classified correctly as NovelInsertion rather than TandemDuplication.

## Success Criteria

- [ ] `cargo test` passes with the fix applied.
- [ ] New test case using benchmark-derived sequences classifies as NovelInsertion.
- [ ] `cargo clippy -- -D warnings` passes.
- [ ] `/update` has been run after changes.

## Steps

1. Apply the fix identified in SV-CLS-001.
2. Add a unit test with real benchmark sequences.
3. Run the full test suite.

## Notes

The fix must not break existing TandemDuplication detection. Ensure the
LargeDeletion and Inversion classification paths are also unaffected.
