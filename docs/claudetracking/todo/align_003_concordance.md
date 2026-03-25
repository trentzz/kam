# ALIGN-003: Per-Variant Concordance Table

**Epic**: ALIGN-COMPARE (docs/claudetracking/overallplans/ALIGN-COMPARE.md)
**Priority**: high
**Depends on**: align_001_parse_thesis.md, align_002_run_titration.md
**Status**: todo

## Goal

Build a concordance table with one row per (variant × sample) pair comparing
kam and alignment pipeline detection. Done looks like:
`docs/benchmarking/align_compare/concordance_table.tsv` with 9,000 rows (375
variants × 24 samples) and columns for detection status from both methods.

## Steps

1. Read `thesis_parsed/alignment_results.csv` and the 24 kam result TSVs.

2. Write `docs/benchmarking/align_compare/scripts/build_concordance.py`:
   - For each (variant × sample) pair:
     - Look up `detected_by_alignment` from alignment results.
     - Look up whether kam called the variant in the corresponding sample TSV
       (PASS filter, matching chrom + pos + ref + alt within ±2 bp).
     - Record: `both_pass`, `kam_only`, `alignment_only`, `both_miss`.
   - Add columns: `kam_filter`, `kam_vaf`, `kam_confidence`, `true_vaf`,
     `alignment_filter`.
   - Write to `concordance_table.tsv`.

3. Compute summary statistics:
   - Per-method sensitivity: fraction of truth variants detected (PASS) across
     all samples.
   - Concordant PASS rate, kam-only rate, alignment-only rate, both-miss rate.
   - Stratified by variant type (SNV, indel, SV) and by VAF bin (0–1%,
     1–10%, >10%).

4. Write summary statistics to `concordance_summary.tsv`.

## Notes

- Join key: `variant_id` (from alignment results) must match the kam result
  TSV variant identifier. Define the join key carefully to avoid mismatches.
- For variants where the alignment pipeline reports a different position than
  truth (off-by-1 due to normalisation), use a ±2 bp tolerance.
- If a sample is missing from either result set (e.g. a failed kam run), mark
  all rows for that sample as `status = missing` and exclude from sensitivity
  calculations.
