# BENCH-SNV-004: Indel VAF Sweep Configs and Data

**Epic**: BENCH-SNV (docs/claudetracking/overallplans/BENCH-SNV.md)
**Priority**: high
**Depends on**: BENCH-SNV-002
**Status**: done

## Goal

Expand the indel synthetic benchmark to cover 0.5%, 1%, 2%, and 5% VAF.
Add truth VCFs and varforge configs for the three missing VAF levels.
Run varforge to generate all FASTQ datasets.

## Steps

1. Add truth_indels_vaf005.vcf, truth_indels_vaf020.vcf, truth_indels_vaf050.vcf
2. Add configs: indel_vaf005.yaml, indel_vaf020.yaml, indel_vaf050.yaml
3. Run varforge for all four VAF levels

## Notes

Same five indels as the existing 1% dataset. Two small deletions (3 bp, 4 bp),
two small insertions (5 bp), one deletion (4 bp). Spread across the 2000 bp
reference at positions 100, 300, 550, 750, 1000.
