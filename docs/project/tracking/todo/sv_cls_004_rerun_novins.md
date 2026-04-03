# SV-CLS-004: Re-run NovelInsertion Benchmarks with Fix

**Epic**: SV-CLASSIFY
**Priority**: high
**Depends on**: SV-CLS-002, SV-BASE-002
**Status**: todo

## Goal

Rebuild kam with the classification fix and re-run all 50 NovelInsertion
datasets to verify correct classification end-to-end.

## Success Criteria

- [ ] All 50 kam_novins results show NovelInsertion (SVTYPE=INS, not TandemDuplication/DUP).
- [ ] Both discovery and tumour-informed outputs are regenerated.
- [ ] `/update` has been run after changes.

## Steps

1. Rebuild kam with the fix from SV-CLS-002.
2. Re-run kam on all 50 NovelInsertion datasets.
3. Verify SVTYPE in all output VCFs.

## Notes

Save the old results for comparison in SV-CLS-005 before overwriting.
