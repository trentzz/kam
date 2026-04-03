# PUB-002: UMI Clustering Benchmark (SRR6794144)

**Epic**: PUB-BENCH (docs/claudetracking/overallplans/PUB-BENCH.md)
**Priority**: high
**Depends on**: pub_001_dataset_download.md, chem_005_non_twist_integration.md
**Status**: todo

## Goal

Run kam on SRR6794144 (simplex 12 bp UMI) in both discovery and tumour-informed
modes. Compare molecule count and UMI clustering behaviour to the published
reference results. Done looks like: a results TSV and a summary note in
`docs/benchmarking/public/umi_clustering/` documenting molecule counts, UMI
collision rate, and any discrepancies from the published benchmark.

## Steps

1. Download the dataset using the script from PUB-001 (if not already done).

2. Create `docs/benchmarking/public/umi_clustering/config.toml`:
   ```toml
   [chemistry]
   preset = "simplex-12bp"
   ```
   Set other parameters to defaults unless the dataset requires adjustment.

3. Build a minimal target FASTA from the published panel regions for this
   sample (or use a k-mer allowlist from known variants if provided).

4. Run kam in discovery mode:
   ```bash
   kam run --config config.toml --r1 SRR6794144_1.fastq.gz \
       --r2 SRR6794144_2.fastq.gz --targets targets.fasta \
       --output results/calls_discovery.vcf
   ```

5. Run kam in tumour-informed mode using known variant positions from the
   benchmark truth set.

6. Extract from the QC JSON:
   - Total read pairs
   - Molecules assembled (after UMI clustering)
   - UMI collision rate (if reported)
   - Molecules per UMI family (distribution)

7. Compare to the published HUMID or UMI-tools results for the same sample.
   Document agreement and discrepancies in
   `docs/benchmarking/public/umi_clustering/RESULTS.md`.

## Notes

- The primary purpose of this benchmark is to validate CHEM-CONFIG (12 bp UMI
  support), not to benchmark variant calling. Focus on molecule assembly
  quality.
- If the sample has no truth VCF, skip tumour-informed mode and report
  discovery results only.
- Record the kam version (git hash) and config in the results directory.
