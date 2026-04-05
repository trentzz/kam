# SV-BASE-001: Complete InvDel Varforge Simulations

**Epic**: SV-BASELINE
**Priority**: critical
**Depends on**: none
**Status**: todo

## Goal

Run varforge for the remaining 48 of 50 InvDel configs to produce simulated
FASTQ pairs and truth VCFs for every VAF level and replicate.

## Success Criteria

- [ ] All 50 sim_invdel_vaf* directories contain R1/R2 FASTQ files and a truth VCF.
- [ ] `/update` has been run after changes.

## Steps

1. Identify which 2 of 50 InvDel configs have already been run.
2. Run varforge for the remaining 48 configs.
3. Verify output completeness across all 50 directories.

## Notes

Configs are in `docs/benchmarking/sv_new/configs/`. Each config produces a
paired FASTQ and truth VCF under `docs/benchmarking/sv_new/results/`.
