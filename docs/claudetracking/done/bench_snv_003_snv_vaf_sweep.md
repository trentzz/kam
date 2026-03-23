# BENCH-SNV-003: SNV VAF Sweep Configs and Data

**Epic**: BENCH-SNV (docs/claudetracking/overallplans/BENCH-SNV.md)
**Priority**: high
**Depends on**: BENCH-SNV-002
**Status**: done

## Goal

Expand the SNV synthetic benchmark to cover 0.5%, 1%, 2%, and 5% VAF.
Add truth VCFs and varforge configs for the three missing VAF levels.
Run varforge to generate all FASTQ datasets.

## Steps

1. Add truth_snvs_vaf005.vcf, truth_snvs_vaf020.vcf, truth_snvs_vaf050.vcf
2. Add configs: snv_vaf005.yaml, snv_vaf020.yaml, snv_vaf050.yaml
3. Update .gitignore for new result directories
4. Run varforge for all four VAF levels

## Notes

Same five SNVs as the existing 1% dataset — only the VAF INFO tag and purity
field change. Purity formula: purity = VAF × 2 (diploid assumption).
