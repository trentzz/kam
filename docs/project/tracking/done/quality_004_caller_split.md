# QUALITY-004: Split kam-call/src/caller.rs into focused modules

**Epic**: QUALITY
**Priority**: low
**Depends on**: none
**Status**: done

## Goal

Split the 2271-line `kam-call/src/caller.rs` into focused modules: `types.rs` (VariantCall, VariantFilter, VariantType, CallSource), `classify.rs` (classify_variant + helpers), `stats.rs` (estimate_vaf, strand_bias_test, compute_confidence), `filter.rs` (assign_filter + CallerConfig).

## Success Criteria

- [x] `caller.rs` is under 500 lines. (mod.rs is 473 lines)
- [x] Each new module has its own `#[cfg(test)] mod tests`. (classify, filter, stats, mod all have tests; types is pure data)
- [x] No public API change: `pub use` re-exports in `lib.rs` preserve existing import paths.
- [x] All tests pass.
- [x] `/update` has been run after changes.

## Steps

1. Extract types into `types.rs`.
2. Extract classification helpers into `classify.rs`.
3. Extract statistical functions into `stats.rs`.
4. Extract filter assignment into `filter.rs`.
5. Verify doc-tests unchanged.

## Notes

Keep `call_variant` as the primary entry point in `caller.rs`.
