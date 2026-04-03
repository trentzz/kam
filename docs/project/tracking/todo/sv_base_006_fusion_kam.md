# SV-BASE-006: Run kam on Fusion Datasets

**Epic**: SV-BASELINE
**Priority**: critical
**Depends on**: SV-BASE-003
**Status**: todo

## Goal

Run kam in both discovery and tumour-informed modes on all 50 Fusion mixed
datasets to produce variant calls.

## Success Criteria

- [ ] All 50 kam_fusion_vaf* directories contain discovery and tumour-informed VCF and TSV output.
- [ ] `/update` has been run after changes.

## Steps

1. Build kam with the latest code.
2. Run kam discovery mode on all 50 Fusion mixed datasets.
3. Run kam tumour-informed mode on all 50 Fusion mixed datasets.
4. Verify output completeness.

## Notes

Fusion detection requires junction k-mers spanning two target regions.
Ensure the target FASTA includes fusion junction sequences.
