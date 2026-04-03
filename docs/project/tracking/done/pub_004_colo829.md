# PUB-004: COLO829 SV Benchmark

**Epic**: PUB-BENCH (docs/claudetracking/overallplans/PUB-BENCH.md)
**Priority**: high
**Depends on**: pub_001_dataset_download.md, sv_exp_009_vcf_new_sv.md
**Status**: todo

## Goal

Run kam on the COLO829 WGS dataset using junction k-mers designed for the 68
validated SVs. Evaluate SV sensitivity per type in tumour-informed mode. Done
looks like: a per-SV sensitivity table in
`docs/benchmarking/public/colo829/RESULTS.md` with one row per validated SV.

## Steps

1. Download the dataset using the script from PUB-001.

2. Download and parse the published COLO829 SV truth set (Lumpy/DELLY
   consensus or the paper supplementary table). Extract the 68 validated SVs:
   - chrom, pos, end, svtype, partner (for BND/fusions)
   - Write to `docs/benchmarking/public/colo829/truth_svs.tsv`.

3. Design junction k-mers for each of the 68 SVs:
   - For each SV, extract a 200 bp window around the breakpoint from the
     reference genome (hg38 or hg19 — use the same build as the dataset).
   - For fusions: extract both breakpoints and generate fusion junction k-mers
     using the utility from SV-EXP-005.
   - Write all k-mer targets to `docs/benchmarking/public/colo829/targets.fasta`.

4. Run kam in tumour-informed mode:
   ```bash
   kam run --r1 ... --r2 ... --targets targets.fasta \
       --target-variants truth_svs_as_vcf.vcf \
       --output results/calls_tumour_informed.vcf
   ```

5. Score against the truth set:
   - True positive: an SV called within ±50 bp of truth breakpoint, matching
     SVTYPE.
   - Compute sensitivity per SVTYPE group (DEL, DUP, INV, BND/Fusion, INS).

6. Document results in `docs/benchmarking/public/colo829/RESULTS.md`.

## Notes

- COLO829 is WGS and likely has no UMIs. As with SEQC2, use a bypass mode or
  document the limitation.
- 50 bp breakpoint tolerance is standard for SV benchmarking on WGS data
  (SURVIVOR default).
- The COLO829 SV truth set varies by publication; use the most cited version
  (Chiang et al. 2015 or the GIAB SV benchmark if available for this cell line).
