# PAPER-008: Supplementary Per-Variant Detailed Tables

**Epic**: PAPER (docs/claudetracking/overallplans/PAPER.md)
**Priority**: high
**Depends on**: align_003_concordance.md, bench_svn_004_run_score.md
**Status**: todo

## Goal

Produce supplementary tables containing per-variant detailed results for every
benchmark in the paper. Done looks like: TSV and formatted PDF/HTML
supplementary tables committed to `docs/paper/supplementary/`, covering
synthetic benchmarks, public datasets, and the alignment concordance analysis.

## Steps

1. Create `docs/paper/supplementary/` directory.

2. Table S1 — Synthetic benchmark full sensitivity table:
   - Source: BENCH-VARFORGE and BENCH-SV-NEW sensitivity TSVs.
   - Columns: variant_class, vaf, rep, discovery_sensitivity, ti_sensitivity.
   - One file per variant class, or one combined file.

3. Table S2 — Public dataset per-variant results:
   - Source: PUB-002, PUB-003, PUB-004 result TSVs.
   - Columns: dataset, variant_id, chrom, pos, ref, alt, svtype, mode, filter,
     detected.
   - One row per variant per dataset.

4. Table S3 — Concordance table (375 variants × 24 samples):
   - Source: `concordance_table.tsv` from ALIGN-003.
   - This is likely too large for a PDF table; produce as a TSV only with a
     summary version (top 50 discordant variants) in PDF.

5. Table S4 — Discordance analysis:
   - Source: `discordance_table.tsv` from ALIGN-004.
   - All systematically discordant variants with root cause category.

6. Write `docs/paper/supplementary/SUPPLEMENTARY.md` as an index listing
   all tables, their sources, and their column definitions.

7. Convert TSVs to formatted HTML or PDF using a simple Pandoc or Python
   script if the journal requires a formatted supplementary document.

## Notes

- All TSVs must use tab separation and have a header row.
- Column definitions must be written in plain language in the index document.
- If any table exceeds 10,000 rows, provide a subsetted version for the
  formatted supplementary and link to the full TSV in the repository.
