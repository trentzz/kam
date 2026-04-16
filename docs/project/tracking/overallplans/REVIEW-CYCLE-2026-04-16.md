# Review Cycle — 2026-04-16

## Purpose

Pre-release review for v0.3.0, covering new features: `--junction-sequences`, Docker image, and ML locus difficulty features (DUST score + repeat fraction).

## Findings summary

Five actionable findings from the Opus review (three HIGH, two MEDIUM):

| ID | Priority | Finding | Status |
|----|----------|---------|--------|
| RC-001 | HIGH | `ml_filter` uses hardcoded 0.5 threshold, ignoring ModelMeta.ml_pass_threshold | Done |
| RC-002 | HIGH | patient-sv-monitoring.md references `variant_class` column — correct name is `variant_type` | Done |
| RC-003 | HIGH | Dockerfile comment and docker-compose.yml reference stale version 0.2.0 | Done |
| RC-004 | MEDIUM | Example output in patient guide shows n_molecules_ref = 112000 for junction call (should be 0) | Done |
| RC-005 | MEDIUM | CLI docs document `--chemistry` for `run` subcommand — correct flag is `--chemistry-override` | Done |

## Known limitations (not blocking release)

- Junction sequence calls produce degenerate ML features (n_simplex_fwd_alt = 0, n_simplex_rev_alt = 0). ML scoring on junction calls will be unreliable. Documented as a limitation.
- `--threads` and `--log-dir` flags are accepted but not implemented. Already documented in CLI reference.
- Cross-chemistry generalisation and alignment-based comparison deferred to future work (acknowledged in paper).
- No end-to-end test where junction-sequence reads produce a call — only the graceful no-call path is tested.

## Tasks created

- `todo/rc_001_ml_filter_threshold.md`
- `todo/rc_002_variant_class_column.md` (merged into docs fix)
- `todo/rc_003_docker_version_comments.md` (merged into docs fix)
- `todo/rc_004_n_molecules_ref_example.md` (merged into docs fix)
- `todo/rc_005_chemistry_override_docs.md` (merged into docs fix)
