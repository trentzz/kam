# PAPER-006: Results Section — Runtime Comparison

**Epic**: PAPER (docs/claudetracking/overallplans/PAPER.md)
**Priority**: high
**Depends on**: runtime_003_runtime_figure.md
**Status**: todo

## Goal

Write the runtime results subsection: wall-clock comparison of kam vs the
alignment-based pipeline. Done looks like: a complete draft (≈200 words)
committed to `docs/paper/sections/results_runtime.md`.

## Steps

1. Read `docs/benchmarking/runtime/RESULTS.md` and `kam_times.tsv` and
   `alignment_times.tsv`.

2. Write `docs/paper/sections/results_runtime.md` covering:

   **Total runtime (≈100 words)**
   - State total wall-clock time for each pipeline on a representative sample.
   - Quote speedup ratio: "kam completed in X seconds vs Y seconds for the
     alignment pipeline, a Z-fold speedup."
   - Note machine spec (CPU cores, RAM, SSD vs spinning disk).

   **Per-stage breakdown (≈100 words)**
   - Which stage dominates in each pipeline.
   - For kam: is it assembly, pathfinding, or calling?
   - For alignment: is it mapping (BWA/HUMID) or variant calling (km)?
   - Note any stage where alignment is faster and explain why.

3. If runtime results are not yet available, leave `[TODO: fill after
   RUNTIME epic]` placeholders.

## Notes

- If kam is slower than alignment on any metric, report it. A paper that only
  reports favourable results is less credible.
- Memory usage (peak RSS) is also interesting if it differs substantially.
  Include a sentence if the difference is >2×.
- This section is short and does not need a dedicated figure; the runtime
  comparison figure from RUNTIME-003 is referenced here but detailed in
  paper_007.
