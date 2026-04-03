# ALIGN-005: Alignment vs kam Comparison Figures

**Epic**: ALIGN-COMPARE (docs/claudetracking/overallplans/ALIGN-COMPARE.md)
**Priority**: high
**Depends on**: align_004_discordance.md
**Status**: todo

## Goal

Produce three figures comparing kam and alignment pipeline performance on the
24 titration samples. Done looks like: three PNG files in
`docs/benchmarking/align_compare/figures/` ready for inclusion in the paper.

## Steps

1. Read the graph style guide at `~/.claude/guides/graph-style.md`.

2. Figure 1 — Sensitivity comparison bars:
   - Grouped bar chart. Groups: SNV, indel, SV. Two bars per group: kam and
     alignment. Y-axis: overall sensitivity (fraction of truth variants with
     PASS status, averaged across all 24 samples).
   - Error bars: standard deviation across samples.
   - Save to `figures/sensitivity_comparison.png`.

3. Figure 2 — Per-variant heatmap:
   - 375 variants × 24 samples. Cell colour: `both_pass` (green),
     `kam_only` (blue), `alignment_only` (orange), `both_miss` (grey).
   - Sort variants by VAF bin and type. Sort samples by tumour fraction.
   - Save to `figures/concordance_heatmap.png`.
   - Note: this figure may be better suited for supplementary if it is too
     dense for the main text.

4. Figure 3 — Venn diagram:
   - Total PASS counts: kam only, alignment only, both, neither.
   - Use a proportional area Venn or an UpSet plot if more than 2 sets.
   - Save to `figures/detection_venn.png`.

5. Write `docs/benchmarking/align_compare/scripts/plot_comparison_figures.py`
   that produces all three figures from `concordance_table.tsv`.

## Notes

- Use Australian English for axis labels and titles.
- The heatmap is likely too large for the main paper figure. Produce a
  compact version (50 selected variants × 24 samples) for the main text and
  the full version for supplementary.
- Both VCF and TSV must exist from ALIGN-002/003 before this task runs.
