# BENCH-VF-004: Run All Datasets Through kam

**Epic**: BENCH-VARFORGE (docs/claudetracking/overallplans/BENCH-VARFORGE.md)
**Priority**: high
**Depends on**: BENCH-VF-001, BENCH-VF-002, BENCH-VF-003
**Status**: todo

## Goal

Run kam on every generated dataset in both discovery and tumour-informed modes.
Produce `calls_discovery.vcf` and `calls_tumour_informed.vcf` for each dataset.

## Steps

1. Write `docs/benchmarking/snvindel/scripts/run_snv_indel_suite.sh` to run kam
   on all SNV and indel datasets (both modes per dataset).
2. Write `docs/benchmarking/sv/scripts/run_sv_suite.sh` to run kam on all SV
   datasets (both modes per dataset).
3. Run both scripts. Log failures.
4. Verify output files exist for all datasets.

## Notes

Output directory structure:
```
results/kam_snv_vaf0050_a/
  calls_discovery.vcf
  calls_tumour_informed.vcf
results/kam_indel_vaf0050_a/
  calls_discovery.vcf
  calls_tumour_informed.vcf
results/kam_sv_vaf0050_a/
  calls_discovery.vcf
  calls_tumour_informed.vcf
```

For tumour-informed mode, pass the corresponding truth VCF as `--target-variants`.
Use the same kam parameters as the existing benchmark scripts for comparability.
Results directories are gitignored; only the run scripts and a summary TSV are
committed.
