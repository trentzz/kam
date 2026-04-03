# SV-BASE-005: Run kam on NovelInsertion Datasets

**Epic**: SV-BASELINE
**Priority**: critical
**Depends on**: SV-BASE-002
**Status**: todo

## Goal

Run kam in both discovery and tumour-informed modes on all 50 NovelInsertion
simulated datasets to produce variant calls.

## Success Criteria

- [ ] All 50 kam_novins_vaf* directories contain discovery and tumour-informed VCF and TSV output.
- [ ] `/update` has been run after changes.

## Steps

1. Build kam with the latest code.
2. Run kam discovery mode on all 50 NovelInsertion datasets.
3. Run kam tumour-informed mode on all 50 NovelInsertion datasets.
4. Verify output completeness.

## Notes

Results may initially show TandemDuplication instead of NovelInsertion due to
the classification bug tracked in SV-CLASSIFY. Record results as-is for
before/after comparison.
