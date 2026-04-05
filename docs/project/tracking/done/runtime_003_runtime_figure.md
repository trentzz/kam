# RUNTIME-003: Runtime Comparison Figure

**Epic**: RUNTIME (docs/claudetracking/overallplans/RUNTIME.md)
**Priority**: medium
**Depends on**: runtime_002_alignment_timing.md
**Status**: todo

## Goal

Produce a bar chart comparing per-stage wall-clock runtime for kam vs the
alignment pipeline. Done looks like: a PNG at
`docs/benchmarking/runtime/figures/runtime_comparison.png` showing both
pipelines side by side with stacked stage bars.

## Steps

1. Read the graph style guide at `~/.claude/guides/graph-style.md`.

2. Write `docs/benchmarking/runtime/scripts/plot_runtime.py`:
   - Load `kam_times.tsv` and `alignment_times.tsv`.
   - Produce a grouped stacked bar chart:
     - X-axis: sample (or read pair count if you want to show scaling).
     - Y-axis: wall-clock time (seconds, log scale if range is large).
     - Stack bars by stage: each pipeline stage is a different colour.
     - Two bar groups per x-tick: kam (left) and alignment (right).
   - Add a second panel or annotation showing total pipeline time with speedup
     ratio.

3. Save to `docs/benchmarking/runtime/figures/runtime_comparison.png`.

4. Add a summary note in
   `docs/benchmarking/runtime/RESULTS.md`:
   - Total time for kam vs alignment on the representative samples.
   - Which stage dominates in each pipeline.
   - Speedup ratio (alignment_total / kam_total).

## Notes

- Use Australian English for all labels and annotations.
- If kam is slower than the alignment pipeline for any sample or stage, report
  this honestly. Do not cherry-pick samples.
- If only a single sample size was measured (no scaling data), produce a single
  side-by-side grouped bar chart rather than a line plot.
