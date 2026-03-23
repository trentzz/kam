# BENCH-SV-002: New SV Types VAF Sweep

**Epic**: BENCH-SV (docs/claudetracking/overallplans/BENCH-SV.md)
**Priority**: medium
**Depends on**: BENCH-SV-001
**Status**: done

## Goal

Expand the INS, large DEL, and INVDEL benchmarks to cover the same four VAF
levels (0.5%, 1%, 2%, 5%) as the existing DEL/DUP/INV benchmark. The current
new-type configs only test 1% VAF.

## Steps

1. Add truth VCFs for INS at 0.5%, 2%, 5% VAF
2. Add truth VCFs for large DEL (600 bp) at 0.5%, 2%, 5% VAF
3. Add truth VCFs for INVDEL at 0.5%, 2%, 5% VAF
4. Add varforge configs for all nine new combinations
5. Update .gitignore for new result directories
6. Run varforge for all nine new datasets

## Notes

Same variant sequences used at all VAF levels — only the VAF INFO tag and
purity field change. This gives nine new datasets (3 types × 3 new VAF levels)
to complement the three existing 1% VAF datasets.

Large DEL configs at all VAF levels (including the pre-existing 1%) panic in
varforge with an out-of-bounds slice error in engine.rs. This is a known varforge
bug, the same class as the vaf0025 issue. The configs and truth VCFs are committed;
the FASTQs cannot be generated until varforge is fixed upstream. INS and INVDEL
simulate cleanly at all VAF levels.
