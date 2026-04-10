# Changelog

All notable changes to this project are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).
Versioning follows [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.2.0] - 2026-04-11

### Added

- `twist-duplex-v1` ML model bundled with the binary for ML-boosted variant
  scoring in tumour-informed mode.
- `single-strand-v1` ML model bundled for single-strand library support.
- Named ML model resolution: pass a model name (e.g. `twist-duplex-v1`)
  instead of a file path.
- ML scoring propagates `ml_prob` and `ml_filter` columns to TSV output.
- Tumour-informed (TI) mode evaluation on 24-sample titration dataset
  (3 ng × 8 VAF levels, 375-variant truth set).
- Per-sample per-target TSV output from batch benchmark scripts.
- Benchmark `07-snvindel-ml-boost-v1`: ML boost evaluation with full report
  and sensitivity-by-VAF figures.

### Fixed

- ONNX model `ZipMap` node removed from both bundled models; ORT 2.0.0-rc.12
  panicked on the `seq(map(int64, tensor(float)))` output type.
- ML feature compilation: `kam-call` now correctly gates `ort`/`ndarray` on
  the `ml` feature flag, restoring `#[cfg(feature = "ml")]` guards.
- Feature propagation: top-level `kam` crate now passes `kam-call/ml` when
  building with `--features ml`.
- Truth VCF coordinate convention: TI mode now uses 1-based positions
  consistent with Rust internal representation.

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
