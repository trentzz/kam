# DEPS: Dependency and Security Advisory Cleanup

**Status**: todo
**Priority**: medium
**Branch**: epic/DEPS

## Goal

Address three active RustSec advisories (bincode 1.3.3 unmaintained, rand 0.8.5 unsound, paste 1.0.15 unmaintained) and set up `cargo audit` in CI so new advisories surface automatically.

## Motivation

`cargo audit` on 2026-04-20 reports three warnings. Vision goal "Built to production quality: tested, documented, reproducible" requires addressing known-unmaintained deps before paper submission. bincode in particular is load-bearing — all inter-stage file formats use it.

## Child tasks

| ID | File | Status |
|----|------|--------|
| DEPS-001 | todo/deps_001_bincode_decision.md | todo |
| DEPS-002 | todo/deps_002_statrs_rand_upgrade.md | todo |
| DEPS-003 | todo/deps_003_cargo_audit_ci.md | todo |

## Scope

- `Cargo.toml`, `Cargo.lock` across all crates
- `docs/project/devmanual/release.md` — document audit-clean release policy

## Out of scope

- Full migration of on-disk formats (tracked under needs-review/002)
