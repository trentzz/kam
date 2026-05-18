# QUALITY-002: Add GitHub Actions CI

**Epic**: QUALITY
**Priority**: medium
**Depends on**: none
**Status**: done

## Goal

Add `.github/workflows/ci.yml` running `cargo fmt --check`, `cargo clippy --all --all-targets -- -D warnings`, and `cargo test --all` on push to main and on pull requests.

## Success Criteria

- [x] `.github/workflows/ci.yml` exists and passes on main.
- [x] All three jobs (fmt, clippy, test) run on Linux stable Rust.
- [x] README has a CI badge.
- [x] `/update` has been run after changes.

## Steps

1. Write `.github/workflows/ci.yml`.
2. Push to main and verify all jobs green.

## Notes

Defer `cargo audit` job to DEPS-003.
