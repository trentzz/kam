# Paper Figures

Figures are listed in priority order. The first figure is the headline result and must be compelling on its own.

---

## Figure 1: Tumour-Informed Sensitivity Comparison

kam versus alignment-based across VAF levels. This is the headline figure.

- Side-by-side or overlay line plots.
- Show all three input amounts: 5 ng, 15 ng, 30 ng.
- x-axis: VAF level. y-axis: sensitivity.
- Error bars where applicable.

---

## Figure 2: Precision Comparison

kam precision 1.0 versus alignment-based false positive rates.

- Bar chart or summary table.
- Highlight the zero false positive result clearly.

---

## Figure 3: SV Detection Heatmap

Sensitivity by SV type and VAF level.

- Rows: SV types (fusions, translocations, large deletions, inversions, duplications).
- Columns: VAF levels.
- Colour encodes sensitivity.
- Demonstrates the breadth of SV support.

---

## Figure 4: Runtime Comparison

Wall-clock time for kam versus the alignment-based pipeline.

- Bar chart on matched datasets.
- Include error bars if multiple runs are available.

---

## Figure 5: Architecture Diagram

The five-crate pipeline showing molecule flow and information preservation.

- Boxes for each crate (kam-core, kam-assemble, kam-index, kam-pathfind, kam-call).
- Arrows showing data flow.
- Annotate where molecule provenance is preserved.

---

## Figure 6: Per-Variant Analysis

Which variants each method detects.

- Scatter or volcano plot.
- Highlight variants detected by kam but not alignment-based, and vice versa.
- Supports the per-variant deep dive discussion.

---

## Figure 7: Cross-Chemistry Results

Sensitivity on non-Twist datasets.

- Demonstrates generalisation claim.
- Bar chart or line plot by chemistry type.

---

## Figure 8: Molecule Evidence Detail

An example showing how molecule-level tracking provides richer evidence than read-level counting.

- Illustrative: one variant, two methods.
- Show molecule count, duplex support, and strand information alongside a raw read count.

---

## Presentation Preferences

- Clean, minimal style. No chart junk.
- Consistent colour palette across all figures.
- Error bars or confidence intervals wherever applicable.
- Vector format (SVG or PDF) for publication.
