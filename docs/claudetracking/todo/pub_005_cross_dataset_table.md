# PUB-005: Cross-Dataset Results Table

**Epic**: PUB-BENCH (docs/claudetracking/overallplans/PUB-BENCH.md)
**Priority**: high
**Depends on**: pub_002_umi_clustering.md, pub_003_seqc2.md, pub_004_colo829.md
**Status**: todo

## Goal

Produce a single cross-dataset results table aggregating sensitivity and
precision from all three public benchmarks, in both discovery and
tumour-informed modes. Done looks like: a TSV at
`docs/benchmarking/public/cross_dataset_table.tsv` and a rendered version
suitable for the paper methods section.

## Steps

1. Define the table schema:
   - Columns: dataset, mode, variant_type, n_variants, sensitivity, precision.
   - One row per (dataset × mode × variant_type) combination.

2. Write `docs/benchmarking/public/scripts/build_cross_dataset_table.py`:
   - Load results TSVs from `umi_clustering/`, `seqc2/`, `colo829/`.
   - Standardise columns across datasets.
   - Write to `cross_dataset_table.tsv`.

3. Run the script and verify the table is complete (no missing rows).

4. Produce a formatted version for the paper:
   - Group by dataset, with sub-rows for each variant type.
   - Highlight cells where sensitivity ≥ 0.90 (tumour-informed) and ≥ 0.70
     (discovery) in the paper figure/table.

5. Commit both the raw TSV and the formatted version.

## Notes

- For datasets without precision (e.g. UMI Clustering Benchmark, no truth VCF
  for all positions), leave precision blank and note the reason.
- If any dataset produced zero detections, include it in the table as
  sensitivity = 0.0 with a note.
- The table will be referenced directly in paper_005 (cross-chemistry results).
