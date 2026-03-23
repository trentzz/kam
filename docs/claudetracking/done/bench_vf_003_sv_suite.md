# BENCH-VF-003: SV Config Suite (50 configs per working type)

**Epic**: BENCH-VARFORGE (docs/claudetracking/overallplans/BENCH-VARFORGE.md)
**Priority**: high
**Depends on**: none
**Status**: todo

## Goal

Generate 50 varforge configs and matching truth VCFs for each working SV type:
DEL+DUP+INV, INS, and INVDEL. 25 VAF levels × 2 replicates each. Large DEL is
excluded (varforge panic — see BENCH-SV-002 notes). Run varforge for all datasets.

## Steps

1. Write `docs/benchmarking/sv/scripts/make_sv_suite.py` to generate truth VCFs
   and configs for DEL+DUP+INV, INS, and INVDEL across all 25 VAF × 2 replicates.
2. Run to create all files in `data/` and `configs/`.
3. Run varforge for all configs (three type groups, can parallelise across types).
4. Verify all datasets generated without error. Note any VAF levels that trigger
   the engine.rs panic (currently only large DEL is known to fail).

## Notes

SV types use the existing truth variant sequences from `truth_svs_vaf010.vcf`,
`truth_ins_vaf010.vcf`, and `truth_invdel_vaf010.vcf`.
Config naming: `sv_vaf0050_a.yaml`, `ins_vaf0050_a.yaml`, `invdel_vaf0050_a.yaml`.
Truth VCF naming follows the same pattern.
Very low VAF levels (< 0.25%) may trigger varforge panics — test 0.05% and
0.075% first and note any failures before running the full suite.
