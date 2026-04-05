# ML-BOOST: Gradient Boosting for Variant Call Filtering

**Status**: active
**Priority**: medium
**Branch**: epic/ML-BOOST

## Goal

Investigate whether LightGBM or XGBoost can improve variant call quality in
kam. The current filter is rule-based (strand bias p-value, confidence score,
NotTargeted). A trained classifier could reduce false positives at low VAF and
provide a more principled call score than the current heuristics.

The deliverable is a proof-of-concept: train a classifier on the existing
varforge benchmark data, measure improvement in precision/recall on held-out
samples, write a report on findings, and if results are positive, design a
concrete `--ml-model` CLI integration.

See `docs/research/ml_gradient_boosting.md` for the full research report.

## Motivation

At very low VAF (0.05%--0.25%), the pipeline makes calls with nalt=1 or
nalt=2 where confidence is driven purely by molecule count. A classifier
trained on the full feature set (vaf, nref, nalt, ndupalt, nsimalt, sbp,
conf, derived features) may learn a better boundary than the current
threshold rules.

## Design

### Data sources

- `docs/benchmarking/snvindel/samples/` — 100 samples (SNV + indel), each
  with `discovery.tsv`, `tumour_informed.tsv`, `truth.tsv`
- `docs/benchmarking/sv/samples/` — 150 samples (SV: DEL/DUP/INV + INS +
  INVDEL), each with the same structure

### Feature set

Raw: vaf, nref, nalt, ndupalt, nsimalt, sbp, conf, ref_len, alt_len
Derived: duplex_frac (ndupalt/nalt), has_duplex (bool), ci_width (vaf_hi -
vaf_lo), variant_class (SNV/indel/SV from allele lengths)

### Labels

Match each call to the truth VCF by chrom+pos+ref+alt. True positive = 1,
false positive = 0. Ignore NotTargeted calls when evaluating tumour-informed
mode.

### Training strategy

Group k-fold CV with sample as the grouping key. Never split within a sample.
Use 5-fold CV (20 samples per fold for SNV/indel). Evaluate: precision,
recall, F1, AUPRC.

### Integration target

A `--ml-model <path>` flag on the `kam call` subcommand. Emits `ml_prob` and
`ml_filter` columns alongside existing output. Does not remove calls, only
annotates.

## Child tasks

| ID | File | Status |
|----|------|--------|
| ML-001 | todo/ml_001_data_pipeline.md | todo |
| ML-002 | todo/ml_002_train_eval.md | todo |
| ML-003 | todo/ml_003_rust_integration_design.md | todo |
| ML-004 | todo/ml_004_figures_report.md | todo |
