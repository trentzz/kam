# BENCH-SVN-005: Sensitivity Figures by New SV Type

**Epic**: BENCH-SV-NEW (docs/claudetracking/overallplans/BENCH-SV-NEW.md)
**Priority**: high
**Depends on**: bench_svn_004_run_score.md
**Status**: todo

## Goal

Produce sensitivity curve figures for InvDel, NovelInsertion, and Fusion,
following the same style as the existing SV sensitivity figures from
BENCH-VARFORGE. Done looks like: three PNG figures exist in
`docs/benchmarking/sv_new/results/figures/`, one per SV type, each showing
discovery and tumour-informed sensitivity across VAF levels.

## Steps

1. Read `docs/claudeguide/benchmarking.md` and the existing figure scripts
   in `docs/benchmarking/sv/scripts/` to match the figure style.

2. Write `docs/benchmarking/sv_new/scripts/plot_sv_new_sensitivity.py`:
   - Load `sensitivity_invdel.tsv`, `sensitivity_novins.tsv`, and
     `sensitivity_fusion.tsv` from `results/`.
   - For each SV type, produce a figure with:
     - x-axis: VAF (%) on log scale.
     - y-axis: sensitivity (0–1).
     - Two lines: discovery mode (dashed) and tumour-informed mode (solid).
     - Scatter points for individual replicates (A and B).
     - Shaded region between replicate A and B to show stochastic variance.
   - Title: e.g. "InvDel Sensitivity vs VAF".

3. Run the script and save figures to
   `docs/benchmarking/sv_new/results/figures/sensitivity_{type}.png`.

4. Verify figures are readable and axes are correctly labelled.

5. Add a combined 3-panel figure (`sensitivity_all_sv_new.png`) with one
   panel per SV type for use in the paper.

## Notes

- Follow the graph style guide at `~/.claude/guides/graph-style.md`.
- Use the same colour scheme and line styles as the existing BENCH-VARFORGE
  figures so results are visually comparable across epics.
- Both VCF and TSV outputs must exist from BENCH-SVN-004 before this task
  runs, as per the memory note on always outputting both formats.
