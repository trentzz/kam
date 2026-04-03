# PAPER-007: Discussion Section and Final Figures

**Epic**: PAPER (docs/claudetracking/overallplans/PAPER.md)
**Priority**: high
**Depends on**: paper_003_results_ti.md, paper_004_results_sv.md, paper_005_results_chemistry.md, paper_006_results_runtime.md
**Status**: todo

## Goal

Write the discussion section and assemble all final paper figures. Done looks
like: `docs/paper/sections/discussion.md` is complete, all five main figures
exist as high-resolution PNGs in `docs/paper/figures/`, and figure captions
are written.

## Steps

1. Write `docs/paper/sections/discussion.md` (≈600 words) covering:

   **Strengths (≈200 words)**
   - Alignment-free: no reference bias, no mapping artefacts.
   - Molecule-level evidence: UMI provenance preserved throughout.
   - Multi-chemistry: works with any UMI length.
   - Tumour-informed mode: near-zero FP with high sensitivity.

   **Limitations (≈200 words)**
   - Targeted mode only: requires a target FASTA.
   - De novo discovery not yet implemented.
   - Fusion detection requires pre-specified partner pairs.
   - k-mer size limits sensitivity for variants in repetitive regions.

   **Future work (≈100 words)**
   - De novo discovery (Phase 2).
   - Nextflow integration for production use.
   - Larger public dataset validation.

   **Conclusion (≈100 words)**
   - One-paragraph summary: what was built, what it can do, why it matters.

2. Finalise all paper figures:
   - Figure 1: pipeline architecture diagram (commission or create in draw.io).
   - Figure 2: sensitivity curves, tumour-informed mode (from BENCH-VARFORGE).
   - Figure 3: SV sensitivity heatmap (from BENCH-SV-NEW).
   - Figure 4: cross-dataset results table as figure (from PUB-BENCH).
   - Figure 5: runtime comparison bar chart (from RUNTIME).

3. Write figure captions in `docs/paper/figures/CAPTIONS.md`.

4. Commit all figures as PNG (≥300 DPI).

## Notes

- Discussion should interpret results, not repeat them. Each paragraph has
  a clear interpretive point.
- Limitations section must be honest and complete. Reviewers find what
  authors omit.
- Write in Australian English, active voice, no em dashes, no padding.
