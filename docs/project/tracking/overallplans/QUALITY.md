# QUALITY: Library-Code Error Handling and Code Health

**Status**: todo
**Priority**: medium
**Branch**: epic/QUALITY

## Goal

Bring kam source in line with CLAUDE.md rule "No `unwrap()` in library code — use proper error handling with `thiserror`". Remove dead code, uncommitted scratch files, and stale TODOs. Set up CI.

## Motivation

322 `unwrap()` calls across 24 source files directly violate the project's stated code quality rules. 18 uncommitted summary TSVs clutter the working tree. No CI means regressions can land unnoticed. These are small, high-leverage fixes for a mature-looking codebase.

## Child tasks

| ID | File | Status |
|----|------|--------|
| QUALITY-001 | todo/quality_001_unwrap_to_thiserror.md | todo |
| QUALITY-002 | todo/quality_002_ci_workflow.md | todo |
| QUALITY-003 | todo/quality_003_uncommitted_cleanup.md | todo |
| QUALITY-004 | todo/quality_004_caller_split.md | todo |
| QUALITY-005 | todo/quality_005_targeting_dead_else_if.md | todo |

## Scope

- All `src/` files, particularly `kam-assemble/src/parser.rs` and `kam/src/commands/run.rs`
- `.github/workflows/` — new CI config
- `docs/benchmarking/01-snvindel/summary/` — triage stale TSVs

## Out of scope

- `kam-ml` vs `kam-call` consolidation (see needs-review/003)
