# BENCH-SVN-003: Fusion Varforge Configs and Truth VCFs

**Epic**: BENCH-SV-NEW (docs/claudetracking/overallplans/BENCH-SV-NEW.md)
**Priority**: high
**Depends on**: sv_exp_005_fusion_junctions.md, sv_exp_006_fusion_classify.md
**Status**: todo

## Goal

Generate 50 varforge configs and matching truth VCFs for Fusion benchmarking:
25 VAF levels × 2 replicates. Done looks like: all 50 configs and truth VCFs
exist, varforge produces FASTQ for all 50, and the junction k-mer FASTA for
the fusion partner pair is committed alongside the configs.

## Steps

1. Design the fusion partner pair to benchmark:
   - Partner A: positions 200–400 on the 2000 bp synthetic reference (region A).
   - Partner B: positions 1200–1400 on the 2000 bp synthetic reference (region B).
   - Fusion breakpoint: end of region A + start of region B.

2. Generate the junction k-mer FASTA using the utility from SV-EXP-005:
   - Write it to `docs/benchmarking/sv_new/data/fusion_junctions.fasta`.
   - Use k = 31 (same as the index default).

3. Write `docs/benchmarking/sv_new/scripts/make_fusion_suite.py`:
   - 25 VAF levels × 2 replicates.
   - Config naming: `fusion_vaf{tag}_{rep}.yaml`.
   - Truth VCF naming: `truth_fusion_vaf{tag}_{rep}.vcf`.
   - Each config references `fusion_junctions.fasta` as the targets file.

4. Run the script to generate all 50 configs and truth VCFs.

5. Run varforge for all 50 configs. Verify all FASTQ pairs are produced.

## Notes

- Varforge must support simulating chimeric reads that span two regions of the
  same reference. If this is not natively supported, construct the chimeric
  reads manually using a Python script that concatenates read prefixes from
  region A and region B.
- If varforge cannot simulate fusions on a single reference, use two separate
  synthetic references (one per partner) and concatenate reads from both.
- The truth VCF for a fusion contains two BND records per event, matching the
  VCF output format from SV-EXP-009.
