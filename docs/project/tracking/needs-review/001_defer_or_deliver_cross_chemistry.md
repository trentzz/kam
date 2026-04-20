# Cross-Chemistry and Alignment Comparison: Deliver Before Submission or Change Paper Framing

**Category**: scope
**Related epic**: PAPER, PUB-BENCH, ALIGN-COMPARE

## Context

The current paper draft (`docs/paper/sections/05_results.tex`) defers two results the vision declares load-bearing:

- `vision/paper/focus.md` lists "The approach generalises across UMI chemistries" as a secondary claim. Paper section 5.6 is a single paragraph: "Generalisation to other UMI chemistries is reserved for future work."
- `vision/results/comparisons.md` lists the thesis RaSCALL alignment baseline as the primary comparison. Paper section 5.5 is also three lines: "a head-to-head comparison is reserved for future work."

This materially weakens two of the seven figures in `vision/paper/figures.md` (Figure 1 and Figure 7). The paper's current headline is "kam achieves precision 1.0 in tumour-informed mode" — true, but the vision wants "comparable to alignment-based". Without the comparison, "comparable" is unsupported.

## Options

1. **Deliver**. Complete PUB-BENCH (at least UMI-benchmark dataset running end-to-end after CHEM-001), and complete ALIGN-COMPARE concordance analysis. Roughly 3 weeks of work gated on CHEM-CONFIG completion.
2. **Reframe the paper**. Drop the "generalises" and "comparable to alignment-based" claims. Paper becomes "kam achieves precision 1.0 in Twist tumour-informed mode, with molecule-level evidence as a differentiator". Honest; narrower.
3. **Partial deliver**. Run UMI-benchmark only (1 public dataset, 1 chemistry) and parse existing thesis RaSCALL numbers into a concordance table without rerunning anything.

## Recommendation

Option 3. The thesis RaSCALL results are in `docs/benchmarking/04-comparison/alignment_baseline.csv` already — a concordance analysis is python-only work (1–3 days). UMI-benchmark is the lightest public dataset (simplex 12 bp UMI, no SV complexity). Together they back the two claims without blocking on SEQC2/COLO829 or full CHEM-CONFIG, and still let the paper carry the vision's headline narrative.
