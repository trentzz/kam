# Changelog

All notable changes to this project are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).
Versioning follows [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.2.0] - 2026-04-14

### Added

- `twist-duplex-v2` ML model: real-data retrained LightGBM classifier on 24
  titration samples. Reduces discovery-mode false positives by 79--93% while
  retaining 97--100% sensitivity at 1--2% VAF. Precision at 15ng 1% VAF
  improves from 0.755 (raw discovery) to 0.969 (vs tumour-informed 1.000).
  AUPRC 0.973 on held-out 30ng test set.
- 16 new ML features (33 → 49 total): sequence context (subst_type,
  trinuc_context, is_cpg, gc_content_ref, homopolymer_run), additional
  VariantCall fields (n_simplex_fwd_alt/rev_alt, n_duplex_ref, n_simplex_ref,
  mean_alt_error_prob, min_variant_specific_duplex,
  mean_variant_specific_molecules), and derived strand/duplex features
  (strand_asymmetry_alt, duplex_vaf, simplex_vaf, duplex_simplex_vaf_delta).
  Old models are unaffected (metadata-driven feature lookup).
- `twist-duplex-v1` ML model bundled with the binary for ML-boosted variant
  scoring in tumour-informed mode.
- `single-strand-v1` ML model bundled for single-strand library support.
- Named ML model resolution: pass a model name (e.g. `twist-duplex-v2`)
  instead of a file path.
- ML scoring propagates `ml_prob` and `ml_filter` columns to TSV output.
- `--ti-rescue` diagnostic flag: probes k-mer evidence for variants not
  detected in TI mode, reporting alt k-mer coverage and molecule counts for
  each targeted position.
- Per-sample per-target TSV output from batch benchmark scripts.
- Benchmark `07-snvindel-ml-boost-v1`: ML boost evaluation across 24
  titration samples (3 ng × 8 VAF), comparing baseline, TI-only, and ML+TI.
  Full report and sensitivity-by-VAF analysis at
  `docs/benchmarking/07-snvindel-ml-boost-v1/`.
- Alignment comparison document: kam versus alignment-based detection on the
  titration dataset, with per-variant detection table and gap analysis.
- Benchmarking resources folder (`docs/benchmarking/resources/`) with
  external reference data including the alignment-based titration results.

### Fixed

- Indel TI matching regression: an extraneous `+1` in all VCF position
  formulas in `extract_variant_key`, `deletion_key`, and `insertion_key`
  (`kam-call/src/targeting.rs`) caused the TI filter to compute positions one
  higher than the truth VCF, rejecting all INDEL calls as `NotTargeted`. The
  `+1` was introduced in a squash commit and has been removed. Indel
  sensitivity at 2% VAF (15 ng) returns from ~23 to ~67 TPs out of 170.
- Python scoring normalisation: the same off-by-one was present in the
  `run_titration_batch.py` scoring loop. Normalisation loop condition reverted
  to match the Rust implementation.
- ONNX model `ZipMap` node removed from both bundled models; ORT 2.0.0-rc.12
  panicked on the `seq(map(int64, tensor(float)))` output type.
- ML feature compilation: `kam-call` now correctly gates `ort`/`ndarray` on
  the `ml` feature flag, restoring `#[cfg(feature = "ml")]` guards.
- Feature propagation: top-level `kam` crate now passes `kam-call/ml` when
  building with `--features ml`.
- DFS budget test: corrected `budget_exceeded_returns_partial_and_flag` in
  `kam-pathfind` to use a graph with valid paths to end.

### Changed

- `run_titration_batch.py` (benchmark 07): added `--ti-rescue` support and
  per-sample output. Unused variable removed.

## [0.1.0] - 2026-04-11

### Added

- Full alignment-free variant detection pipeline across five crates:
  `kam-core`, `kam-assemble`, `kam-index`, `kam-pathfind`, `kam-call`, and
  the top-level `kam` binary.
- Molecule-level UMI assembly for Twist duplex chemistry (`5M2S+T` read
  structure, canonical duplex pairing).
- k-mer indexing with molecule provenance, replacing Jellyfish.
- de Bruijn graph pathfinding with DFS expansion, replacing km/kmtools.
- Statistical variant calling with binomial confidence intervals.
- SV detection: insertions, large deletions, and inversion-deletions.
- Tumour-informed (TI) monitoring mode for near-zero false-positive calling.
- ML-based variant scoring via the `kam-ml` crate; bundled models
  `twist-duplex-v1` and `single-strand-v1` ship with the binary.
- Named model resolution so callers specify a model by name rather than path.
- QC flags at the pathfinding stage for downstream filtering.
- Structured QC JSON output at each pipeline stage.
- Pathfind QC flag propagated through calling.

### Fixed

- INS pathfinding hang eliminated via early-exit and predecessor caching.
- DFS expansion budget corrected to prevent runaway graph traversal.
- Molecule counts in pathfinding output now reflect real molecule evidence.
- ML3 truth-matching and training pipeline chromosome format corrected.
- `export_model.py` now loads a pre-trained v3 model instead of retraining.

### Changed

- ML model artefacts consolidated: mislabelled v2 artefacts renamed to `v2b`;
  v3 label reserved for the ML3 training dataset.
- Training CSVs compressed before Nextcloud upload to avoid gateway timeouts.
