# Changelog

All notable changes to this project are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).
Versioning follows [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
