# SV-ULV-002: Run Ultra-Low VAF Varforge Simulations

**Epic**: SV-ULTRAVAF
**Priority**: medium-low
**Depends on**: SV-ULV-001
**Status**: todo

## Goal

Run varforge for all 24 ultra-low VAF configs to produce simulated FASTQ
pairs and truth VCFs.

## Success Criteria

- [ ] All 24 simulation directories contain R1/R2 FASTQ files and truth VCF.
- [ ] `/update` has been run after changes.

## Steps

1. Run varforge for all 24 configs.
2. Verify output completeness.

## Notes

Fusion configs may require the same tumour+WT mixing approach as in
SV-BASE-003. Check whether the mixing script handles sub-0.05% VAF correctly.
