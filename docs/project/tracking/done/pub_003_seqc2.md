# PUB-003: SEQC2 HCC1395 Benchmark

**Epic**: PUB-BENCH (docs/claudetracking/overallplans/PUB-BENCH.md)
**Priority**: high
**Depends on**: pub_001_dataset_download.md
**Status**: todo

## Goal

Run kam on the SEQC2 HCC1395 dataset and evaluate SNV/indel sensitivity and
precision against the SEQC2 Tier 1 truth VCF. Done looks like: per-variant
sensitivity and precision numbers in `docs/benchmarking/public/seqc2/RESULTS.md`,
covering both discovery and tumour-informed modes.

## Steps

1. Download the dataset and truth VCF using the script from PUB-001.

2. Determine whether the dataset has UMIs:
   - If yes: create a `config.toml` with the correct chemistry.
   - If no: investigate whether a bypass mode (`--no-umi` or degenerate single-
     read-per-molecule assembly) is feasible. If not, document the blocker
     and flag for architecture discussion.

3. Design a target FASTA from the SEQC2 Tier 1 truth VCF:
   - Extract unique variant positions.
   - For each position, write a 200 bp window as a target sequence.
   - Write to `docs/benchmarking/public/seqc2/targets.fasta`.

4. Run kam in discovery mode and tumour-informed mode.

5. Score against the Tier 1 truth VCF:
   - True positive: called variant within ±2 bp of truth position, matching
     ref and alt alleles.
   - False negative: truth variant not called.
   - False positive: called variant not in truth VCF.
   - Compute sensitivity and precision.

6. Document in `docs/benchmarking/public/seqc2/RESULTS.md`:
   - Dataset description and truth VCF source.
   - kam config used.
   - Sensitivity and precision per mode (discovery / tumour-informed).
   - Any notable failure modes.

## Notes

- The SEQC2 dataset may not have UMIs. If it does not, this benchmark tests
  the downstream variant calling logic only, not the UMI assembly. Document
  this caveat clearly.
- Use Tier 1 variants only (high-confidence calls with ≥50× depth in both
  tumour and normal). Tier 2 and 3 are excluded.
- Do not attempt to match the SEQC2 WGS results; this is an applications note,
  not a comprehensive benchmarking study.
