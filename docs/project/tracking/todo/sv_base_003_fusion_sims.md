# SV-BASE-003: Complete Fusion Varforge Simulations

**Epic**: SV-BASELINE
**Priority**: critical
**Depends on**: none
**Status**: todo

## Goal

Run all 50 fusion and 50 wild-type varforge simulations, then concatenate each
pair into mixed FASTQs representing the tumour+normal mixture at the target VAF.

## Success Criteria

- [ ] All 50 sim_fusion_vaf*_mixed directories contain sample_R1/R2.fastq.gz.
- [ ] `/update` has been run after changes.

## Steps

1. Run varforge for all 50 fusion configs (tumour component).
2. Run varforge for all 50 matching wild-type configs (normal component).
3. Concatenate fusion + WT FASTQs at the correct ratio for each VAF level.
4. Verify output completeness.

## Notes

Fusion simulation requires mixing tumour and WT reads because varforge does not
natively support translocation-like events at a target VAF.
