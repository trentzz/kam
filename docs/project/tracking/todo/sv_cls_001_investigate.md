# SV-CLS-001: Investigate NovelInsertion Misclassification

**Epic**: SV-CLASSIFY
**Priority**: high
**Depends on**: none
**Status**: todo

## Goal

Extract ref/alt sequences from existing NovelInsertion benchmark results and
trace through `is_tandem_duplication()` to identify why it returns true for
genuine novel insertions.

## Success Criteria

- [ ] Investigation document in `docs/research/` with root cause identified.
- [ ] Document includes example ref/alt sequences and the exact code path that misclassifies.
- [ ] `/update` has been run after changes.

## Steps

1. Extract ref and alt sequences from a kam_novins result VCF.
2. Read `is_tandem_duplication()` in `kam-call/src/caller.rs`.
3. Manually trace the logic with the extracted sequences.
4. Identify the root cause and document it.

## Notes

The function likely checks for substring matches that are too permissive,
causing any insertion longer than the reference to be flagged as a tandem dup.
