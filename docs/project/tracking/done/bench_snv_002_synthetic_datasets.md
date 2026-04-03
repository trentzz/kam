# BENCH-SNV-002: Synthetic SNV and Indel Datasets

**Epic**: BENCH-SNV (`overallplans/BENCH-SNV.md`)
**Priority**: high
**Depends on**: none
**Status**: done

## Goal

Create self-contained varforge configs for synthetic SNV-only and indel-only
benchmarks. Configs and truth VCFs are committed; generated FASTQs are not.

Done when:
- `docs/benchmarking/snvindel/data/` has a small reference FASTA and truth VCFs
  for SNV-only and indel-only variants.
- `docs/benchmarking/snvindel/configs/` has varforge configs for each condition
  (one SNV-only, one indel-only, at 1% VAF).
- `.gitignore` excludes the FASTQ output directories.

## Steps

1. Create `docs/benchmarking/snvindel/data/ref.fa` — synthetic 2000 bp
   reference (similar to the SV benchmark).
2. Create `docs/benchmarking/snvindel/data/truth_snvs_vaf010.vcf` — 5 SNVs on
   the reference at 1% VAF.
3. Create `docs/benchmarking/snvindel/data/truth_indels_vaf010.vcf` — 5 small
   indels (2–10 bp) on the reference at 1% VAF.
4. Create `docs/benchmarking/snvindel/configs/snv_vaf010.yaml` — varforge config
   for SNV-only simulation.
5. Create `docs/benchmarking/snvindel/configs/indel_vaf010.yaml` — varforge
   config for indel-only simulation.
6. Add result directories to `.gitignore`.

## Notes

- Keep the reference short (2000 bp) so varforge runs in seconds.
- Use the same UMI/fragment model as the SV benchmark for consistency.
- The truth VCFs must be self-consistent with the synthetic reference (no
  GRCh38 positions).
