# BENCH-VF-001: SNV Config Suite (50 configs)

**Epic**: BENCH-VARFORGE (docs/claudetracking/overallplans/BENCH-VARFORGE.md)
**Priority**: high
**Depends on**: none
**Status**: todo

## Goal

Generate 50 varforge configs and matching truth VCFs for SNV benchmarking:
25 VAF levels × 2 replicates. All datasets use the existing 2000 bp synthetic
reference and the same 5 SNV positions. Run varforge to produce all FASTQs.

## Steps

1. Write a generation script (`docs/benchmarking/snvindel/scripts/make_snv_suite.py`)
   that produces truth VCFs and YAML configs for all 25 VAF × 2 replicate combinations.
2. Run the script to create all files in `data/` and `configs/`.
3. Run varforge for all 50 configs (can parallelise).
4. Verify all 50 datasets generated without error.

## Notes

VAF levels (see epic for full list): 0.05% through 10.00%.
Replicate A seeds: base_seed + 0. Replicate B seeds: base_seed + 1000.
Base seed formula: `int(vaf * 10000)` (e.g. VAF 0.5% → 50).
Config naming: `snv_vaf0050_a.yaml`, `snv_vaf0050_b.yaml` (four-digit suffix,
e.g. 0050 = 0.50%, 0100 = 1.00%, 1000 = 10.00%).
Truth VCF naming: `truth_snvs_vaf0050_a.vcf`, `truth_snvs_vaf0050_b.vcf`.
