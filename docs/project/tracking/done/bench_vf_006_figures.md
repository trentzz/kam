# BENCH-VF-006: Sensitivity Curve Figures

**Epic**: BENCH-VARFORGE (docs/claudetracking/overallplans/BENCH-VARFORGE.md)
**Priority**: medium
**Depends on**: BENCH-VF-005
**Status**: todo

## Goal

Generate sensitivity curve figures from the aggregate results. One figure per
variant class showing discovery vs tumour-informed sensitivity across the full
VAF range, with replicate scatter to show variance.

## Steps

1. Write `docs/benchmarking/snvindel/scripts/plot_snv_indel_suite.py` to produce:
   - SNV sensitivity curve (discovery vs tumour-informed, with replicate scatter)
   - Indel sensitivity curve (same format)
   - Combined SNV + indel panel (two subplots, single-column width)
2. Write `docs/benchmarking/sv/scripts/plot_sv_suite.py` to produce:
   - Per-SV-type sensitivity curves (DEL+DUP+INV, INS, INVDEL)
3. Export all figures as PDF + PNG at 600 DPI. Single-column width (84 mm).
4. Commit all figure files to `docs/benchmarking/*/results/graphs/`.

## Notes

Follow the graph style guide (`~/.claude/guides/graph-style.md`):
- Okabe-Ito palette
- Discovery: solid line; tumour-informed: dashed line
- Replicate scatter: small transparent points
- x-axis: VAF (%), log scale from 0.05% to 10%
- y-axis: sensitivity (%), 0–100%
- Axis labels with units on every figure
One takeaway per figure — state it in the caption.
