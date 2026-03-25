# ALIGN-001: Parse Thesis Results into Standardised Per-Variant CSV

**Epic**: ALIGN-COMPARE (docs/claudetracking/overallplans/ALIGN-COMPARE.md)
**Priority**: high
**Depends on**: none
**Status**: todo

## Goal

Parse the thesis alignment pipeline result tables (Excel or TSV) into a
standardised per-variant CSV with one row per (variant × sample) pair. Done
looks like: `docs/benchmarking/align_compare/thesis_parsed/alignment_results.csv`
exists with columns for variant ID, sample, VAF, detection status, and filter.

## Steps

1. Locate the thesis result files. See `~/.claude/projects/.../memory/reference_thesis.md`
   for the PDF location; the accompanying data files (Excel/TSV) should be in
   the same directory or referenced in the notes file.

2. Write `docs/benchmarking/align_compare/scripts/parse_thesis_results.py`:
   - Read the Excel/TSV file(s).
   - For each variant × sample row, extract:
     - `variant_id`: unique identifier (chrom_pos_ref_alt or paper ID).
     - `sample`: titration sample name.
     - `true_vaf`: expected VAF for this sample.
     - `detected_by_alignment`: bool (PASS or equivalent).
     - `alignment_filter`: filter value if available.
   - Write to `alignment_results.csv`.

3. Validate the output:
   - Expected row count: 375 variants × 24 samples = 9,000 rows (or subset
     if not all combinations are in the thesis).
   - Check for missing values and document any.

4. Write a brief data description in
   `docs/benchmarking/align_compare/thesis_parsed/README.md`:
   - Source file(s) used.
   - Which table/sheet in the thesis data.
   - Any transformations applied (e.g. filter remapping).

## Notes

- This task can start immediately with no code dependencies.
- If the thesis data is in PDF tables only (not machine-readable), use a PDF
  extraction tool or manual transcription. Document the extraction method.
- `variant_id` must be a stable, unique key that can be joined against kam
  results in ALIGN-003.
