# ALIGN-002: Run kam on All 24 Titration Samples

**Epic**: ALIGN-COMPARE (docs/claudetracking/overallplans/ALIGN-COMPARE.md)
**Priority**: high
**Depends on**: none
**Status**: todo

## Goal

Run `kam run` in tumour-informed mode on all 24 titration samples from the
thesis, producing per-variant results for each sample. Done looks like:
24 result directories in `docs/benchmarking/align_compare/kam_results/`, each
containing `calls_tumour_informed.vcf` and `calls_tumour_informed.tsv`.

## Steps

1. Locate the 24 titration sample FASTQs. The samples are from the Twist panel;
   see `docs/research/twist_umi_chemistry.md` for chemistry details.

2. Write `docs/benchmarking/align_compare/scripts/run_titration_samples.sh`:
   - Array of sample names and corresponding R1/R2 paths.
   - For each sample, run:
     ```bash
     kam run --r1 $R1 --r2 $R2 --targets targets.fasta \
         --target-variants truth.vcf \
         --output kam_results/${sample}/calls_tumour_informed.vcf \
         --tsv kam_results/${sample}/calls_tumour_informed.tsv
     ```
   - Use the Twist UMI duplex config (default settings).
   - Log output to `kam_results/${sample}/kam.log`.

3. Run the script (preferably on HPC with parallelism if available).

4. Verify all 24 result directories are populated and VCFs are non-empty.

5. Record the kam version (git hash) in
   `docs/benchmarking/align_compare/kam_results/RUN_INFO.md`.

## Notes

- Both VCF and TSV outputs are required (per the project memory note on output
  formats).
- Use the tumour-informed mode only for this comparison (matching the thesis
  pipeline's tumour-informed filtering approach).
- If any sample fails, document the error in `RUN_INFO.md` and exclude it from
  the concordance analysis with a note.
