# BENCH-SNV-005: Combined SNV + Indel Datasets

**Epic**: BENCH-SNV (docs/claudetracking/overallplans/BENCH-SNV.md)
**Priority**: medium
**Depends on**: BENCH-SNV-003, BENCH-SNV-004
**Status**: done

## Goal

Create combined SNV + indel truth VCFs and varforge configs at 0.5%, 1%, 2%,
and 5% VAF. A combined dataset is more realistic than single-type datasets and
tests that SNV and indel scoring do not interfere.

## Steps

1. Add truth_combined_vaf{005,010,020,050}.vcf (merge SNV + indel variants)
2. Add configs: combined_vaf{005,010,020,050}.yaml
3. Run varforge for all four VAF levels

## Notes

The combined VCF merges the five SNVs and five indels from the per-type VCFs,
sorted by position. Ten variants total per dataset. Tests the full SNV/indel
calling pathway in a single run.
