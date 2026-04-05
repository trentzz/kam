# ALIGN-COMPARE: Comparison Against Thesis Alignment-Based Results

**Status**: todo
**Priority**: high
**Branch**: epic/ALIGN-COMPARE

## Goal

Produce a systematic per-variant concordance analysis between kam's tumour-informed
mode and the alignment-based RaSCALL results from the thesis. For all 375
variants across 24 titration samples, determine which method detects each
variant, which misses it, and why. Produce summary figures suitable for a
methods paper.

## Motivation

The thesis established a near-zero false-positive rate for the alignment
pipeline using tumour-informed filtering. kam must be compared against this
baseline. A per-variant concordance table is more informative than aggregate
sensitivity numbers: it shows whether kam's misses are random (stochastic at
low VAF) or systematic (a class of variants kam cannot detect).

## Design

### Data sources

- Thesis result tables: per-variant detection status in all 24 samples.
  Location: thesis Excel/TSV files (parse in pub_001 / align_001).
- kam results: run on all 24 titration samples in tumour-informed mode
  (align_002).

### Concordance table

One row per (variant × sample) pair. Columns:
- variant ID, sample, true VAF, expected depth
- detected_by_kam (bool), detected_by_alignment (bool)
- kam_filter, kam_vaf, kam_confidence
- alignment_filter (PASS/LowDP/etc.)

From this table, compute:
- Per-variant sensitivity for each method
- Concordant PASS, kam-only PASS, alignment-only PASS, both-miss rates
- Stratified by variant type (SNV, indel, SV)

### Discordance analysis

For each discordant case, categorise the root cause:
- VAF below detection threshold for one method
- Alignment artefact (soft-clip, mapping quality)
- k-mer graph artefact (anchor failure, path not found)
- Systematic (affects all samples for this variant)

### Figures

1. Side-by-side sensitivity bars: kam vs alignment, per variant type.
2. Per-variant heatmap: 375 variants × 24 samples, coloured by detection status.
3. Venn diagram: PASS counts per method.

## Child tasks

| ID | File | Status |
|----|------|--------|
| ALIGN-001 | todo/align_001_parse_thesis.md | todo |
| ALIGN-002 | todo/align_002_run_titration.md | todo |
| ALIGN-003 | todo/align_003_concordance.md | todo |
| ALIGN-004 | todo/align_004_discordance.md | todo |
| ALIGN-005 | todo/align_005_comparison_figures.md | todo |

## Dependencies

- ALIGN-001 can start immediately (parse thesis data, no kam changes needed).
- ALIGN-002 requires access to the 24 titration sample FASTQs.
- ALIGN-003 depends on ALIGN-001 and ALIGN-002.

## Scope

- `docs/benchmarking/align_compare/` — scripts, concordance tables, figures
- `docs/benchmarking/align_compare/thesis_parsed/` — standardised thesis TSVs
- `docs/benchmarking/align_compare/kam_results/` — per-sample TSVs from kam
- `docs/benchmarking/align_compare/figures/` — output figures

## Out of scope

- Re-running the alignment pipeline from raw data
- Changes to the kam Rust code
- Statistical models beyond sensitivity/precision counts
