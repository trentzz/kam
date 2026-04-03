# BENCH-VF-002: Indel Config Suite (50 configs)

**Epic**: BENCH-VARFORGE (docs/claudetracking/overallplans/BENCH-VARFORGE.md)
**Priority**: high
**Depends on**: none
**Status**: todo

## Goal

Generate 50 varforge configs and matching truth VCFs for indel benchmarking:
25 VAF levels × 2 replicates. Uses the same 5 indel variants as the existing
1% benchmark. Run varforge to produce all FASTQs.

## Steps

1. Write or extend the generation script to produce indel truth VCFs and configs
   for all 25 VAF × 2 replicate combinations.
2. Run to create all files in `data/` and `configs/`.
3. Run varforge for all 50 configs.
4. Verify all 50 datasets generated without error.

## Notes

Same indels as `truth_indels_vaf010.vcf`: two 3–4 bp deletions, two 5 bp
insertions, one 4 bp deletion. Positions 100, 300, 550, 750, 1000.
Same naming convention as BENCH-VF-001: `indel_vaf0050_a.yaml`, etc.
Can parallelise with BENCH-VF-001 (independent files and directories).
