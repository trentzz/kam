# BENCH-SVN-001: InvDel Varforge Configs and Truth VCFs

**Epic**: BENCH-SV-NEW (docs/claudetracking/overallplans/BENCH-SV-NEW.md)
**Priority**: high
**Depends on**: sv_exp_002_invdel_classification.md
**Status**: todo

## Goal

Generate 50 varforge configs and matching truth VCFs for InvDel benchmarking:
25 VAF levels × 2 replicates. Done looks like: all 50 configs and truth VCFs
exist in `docs/benchmarking/sv_new/configs/` and `docs/benchmarking/sv_new/data/`,
and varforge produces FASTQ for all 50 without error.

## Steps

1. Read the existing SV suite generation script
   (`docs/benchmarking/sv/scripts/make_sv_suite.py` or equivalent) to
   understand the pattern for config generation.

2. Write `docs/benchmarking/sv_new/scripts/make_invdel_suite.py`:
   - 25 VAF levels (same as BENCH-VARFORGE).
   - 2 replicates per level (seeds A and B; formula: `int(vaf * 10000)` and
     `int(vaf * 10000) + 1000`).
   - Config naming: `invdel_vaf{tag}_{rep}.yaml`.
   - Truth VCF naming: `truth_invdel_vaf{tag}_{rep}.vcf`.
   - SV definition: 1 InvDel event on the 2000 bp synthetic reference
     (define position and size to match what the classifier can detect).

3. Run the script to generate all configs and truth VCFs.

4. Run varforge for all 50 configs:
   ```bash
   for cfg in docs/benchmarking/sv_new/configs/invdel_*.yaml; do
       varforge run "$cfg"
   done
   ```

5. Verify all 50 FASTQ pairs were generated without error.

## Notes

- The InvDel must be represented in varforge as a deletion with a flanking
  inversion. If varforge does not natively support InvDel, compose it as a
  deletion variant that also flips an adjacent segment.
- Choose an InvDel size (e.g. 80 bp total: 20 bp inverted + 30 bp deleted
  each side) that the classifier can detect reliably at ≥1% VAF.
- The synthetic reference must be the same 2000 bp reference used in other
  SV benchmarks for consistency.
