# PAPER-004: Results Section — SV Detection

**Epic**: PAPER (docs/claudetracking/overallplans/PAPER.md)
**Priority**: high
**Depends on**: bench_svn_005_sv_figures.md
**Status**: todo

## Goal

Write the SV detection results subsection: sensitivity by SV type and VAF
level, covering all five SV types (DEL, DUP, INV, InvDel, NovelInsertion)
and Fusion if benchmarks are complete. Done looks like: a complete draft
(≈400 words) committed to `docs/paper/sections/results_sv.md`.

## Steps

1. Read the SV sensitivity tables from BENCH-VARFORGE and BENCH-SV-NEW.

2. Write `docs/paper/sections/results_sv.md` covering:

   **SV sensitivity overview (≈150 words)**
   - Summary table or prose: sensitivity at 0.5%, 1%, 2% VAF for each SV type.
   - Distinguish discovery vs tumour-informed mode.
   - Note which types achieve >90% sensitivity at what VAF threshold.

   **InvDel and NovelInsertion (≈100 words)**
   - New types added in SV-EXPAND.
   - Sensitivity comparison to existing SV types.

   **Fusion (≈100 words, or note if deferred)**
   - If fusion benchmarks are complete: report sensitivity.
   - If not yet complete: note as planned work with a `[TODO]` placeholder.

   **SV scoring model (≈50 words)**
   - Brief note on `mean_variant_specific_molecules` (from SV-SCORE epic) and
     why it improves sensitivity for large SVs over the minimum-molecule
     bottleneck.

3. Leave `[TODO: insert numbers]` where results are pending.

## Notes

- SV sensitivity is the secondary headline result. Frame it as extending the
  alignment-free approach to structural variants, which alignment pipelines
  struggle with at low VAF.
- Do not explain the scoring model in depth here; the methods section covers
  it. One sentence is enough.
