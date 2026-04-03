# FIX-002: Replace DefaultHasher with a deterministic hasher

**Epic**: none (standalone fix)
**Priority**: low
**Depends on**: none
**Status**: done

## Goal

Replace `std::collections::hash_map::DefaultHasher` in
`kam-assemble/src/assembler.rs` with a hasher that is guaranteed to produce
stable output across Rust versions. Done looks like: molecule IDs derived from
UMI pair hashes are reproducible regardless of compiler version.

## Steps

1. Read `kam-assemble/src/assembler.rs`, finding `hash_umi_pair` (around
   line 492).

2. Replace `DefaultHasher` with a portable alternative. Options:
   - Use `std::hash::SipHasher` with fixed keys (deprecated but stable output).
   - Add `fnv` crate (lightweight, no_std compatible) and use `FnvHasher`.
   - Use a manual FNV-1a implementation (6 lines, no dependency).
   The FNV-1a approach is preferred: no new dependency, deterministic, fast.

3. Update the `Cargo.toml` if a new dependency is added.

4. Add a test that the hash of a fixed input is a fixed value across runs
   (regression test for stability).

5. Move this task file to done and update Status.

6. Run `cargo test` and `cargo clippy -- -D warnings`. Fix all failures.

## Notes

- From codebase review finding NEW-007.
- Project rule: "All RNG must use seeded deterministic generators for Nextflow
  cache compatibility." Hashes used for molecule IDs fall under this rule.
