# PAPER-003: Results Section — Tumour-Informed Performance

**Epic**: PAPER (docs/claudetracking/overallplans/PAPER.md)
**Priority**: high
**Depends on**: align_003_concordance.md
**Status**: todo

## Goal

Write the headline results subsection: tumour-informed mode sensitivity and
precision on the synthetic VAF sweep and the 24 titration samples. Done looks
like: a complete draft (≈500 words) committed to
`docs/paper/sections/results_ti.md`, with concrete numbers from the benchmarks.

## Steps

1. Read the aggregated sensitivity tables from BENCH-VARFORGE
   (`docs/benchmarking/snvindel/results/`) and the concordance summary from
   ALIGN-003 (`docs/benchmarking/align_compare/concordance_summary.tsv`).

2. Write `docs/paper/sections/results_ti.md` covering:

   **SNV/indel sensitivity curve (≈150 words)**
   - Sensitivity at key VAF levels in tumour-informed mode.
   - Quote: "At 0.5% VAF, sensitivity was X% for SNVs and Y% for indels."
   - State the VAF at which sensitivity drops below 90% (the detection
     threshold).

   **False positive rate in tumour-informed mode (≈100 words)**
   - Number of false positives per sample in tumour-informed mode.
   - Compare to discovery mode FP rate to show the value of tumour-informed
     filtering.

   **Titration sample concordance (≈200 words)**
   - Overall concordance with alignment pipeline.
   - Variants detected by both, detected by kam only, detected by alignment
     only.
   - Key finding: whether kam achieves comparable sensitivity to alignment
     in tumour-informed mode.

3. Leave `[TODO: insert numbers]` placeholders where benchmark results are
   not yet available.

## Notes

- This is the headline result. Write it to highlight the core claim: kam
  achieves near-alignment sensitivity in tumour-informed mode without alignment.
- Use active voice: "kam detected 94% of SNVs at 0.5% VAF", not "94% of SNVs
  were detected".
- Do not include figures in this file; figures are referenced by name
  (e.g. "Figure 2") and produced in paper_007.
