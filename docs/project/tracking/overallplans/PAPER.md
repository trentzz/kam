# PAPER: Proof-of-Concept Methods Paper

**Status**: done
**Priority**: high
**Branch**: epic/PAPER

## Goal

Write and submit a proof-of-concept methods paper describing kam: its
architecture, benchmarking results, and comparison to the alignment-based
baseline. The paper covers methods, results, discussion, figures, and
supplementary material. Target venue: Bioinformatics or Briefings in
Bioinformatics (short applications note, ~3000 words + figures).

## Motivation

The thesis established the clinical context and alignment-based baseline. kam
is the alignment-free replacement. A standalone paper separates the tool
contribution from the thesis and makes it citable and independently evaluable.

## Design

### Sections

| Section | Task | Core content |
|---------|------|-------------|
| Methods: tool | paper_001 | Architecture, config, chemistry |
| Methods: benchmarking | paper_002 | Datasets, metrics, comparison design |
| Results: tumour-informed | paper_003 | Sensitivity and precision, headline |
| Results: SV detection | paper_004 | Sensitivity by SV type and VAF |
| Results: cross-chemistry | paper_005 | Public dataset generalisation |
| Results: runtime | paper_006 | Speed comparison |
| Discussion + figures | paper_007 | Interpretation, limitations, figures |
| Supplementary | paper_008 | Per-variant detailed tables |

### Figures (target set)

1. Pipeline architecture diagram (methods).
2. Sensitivity curve: tumour-informed mode across VAF sweep (SNV + indel).
3. SV sensitivity heatmap: type × VAF.
4. Cross-dataset results table (rendered as figure).
5. Runtime bar chart.
6. Supplementary: per-variant concordance heatmap (375 × 24).

### Style

Target journal: Bioinformatics applications note. Max 3000 words main text,
5 figures, unlimited supplementary. Australian English throughout.

## Child tasks

| ID | File | Status |
|----|------|--------|
| PAPER-001 | done/paper_001_methods_tool.md | done |
| PAPER-002 | done/paper_002_methods_bench.md | done |
| PAPER-003 | done/paper_003_results_ti.md | done |
| PAPER-004 | done/paper_004_results_sv.md | done |
| PAPER-005 | done/paper_005_results_chemistry.md | done |
| PAPER-006 | done/paper_006_results_runtime.md | done |
| PAPER-007 | done/paper_007_discussion_figures.md | done |
| PAPER-008 | done/paper_008_supplementary.md | done |

## Dependencies

All results epics must be complete before writing results sections:
- BENCH-VARFORGE (done)
- BENCH-SV-NEW
- PUB-BENCH
- ALIGN-COMPARE
- RUNTIME

Methods sections (paper_001, paper_002) can be drafted in parallel with
benchmarking.

## Scope

- `docs/paper/` — manuscript LaTeX or Markdown source
- `docs/paper/figures/` — final figure files
- `docs/paper/supplementary/` — supplementary tables

## Out of scope

- Journal submission mechanics
- Peer review responses
- Data deposition (Zenodo/SRA)
