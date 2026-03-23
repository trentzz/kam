# FIX-004: classify_variant returns Snv for identical ref and alt sequences

**Epic**: none (standalone fix)
**Priority**: low
**Depends on**: none
**Status**: done

## Goal

Prevent `classify_variant` from returning `VariantType::Snv` when ref and alt
sequences are identical (diffs == 0). Done looks like: identical sequences
return a dedicated `NoVariant` variant (or are handled upstream), and no
spurious 100% VAF calls appear in output for no-op paths.

## Steps

1. Read `kam-call/src/caller.rs`, focusing on `classify_variant` and the
   `diffs == 0` arm (around line 439).

2. Options:
   a. Add `VariantType::NoVariant` and return it when sequences are identical.
      Update all `match` arms on `VariantType` accordingly.
   b. In `call_variant`, skip calling `classify_variant` when sequences are
      identical (short-circuit before the call and return a filtered result).
   Option (b) is simpler and avoids touching the `VariantType` enum (which
   is a core type). Prefer (b) unless `NoVariant` is needed elsewhere.

3. Add a test: `call_variant` with identical ref and alt should produce a
   filtered call (not PASS).

4. Move this task file to done and update Status.

5. Run `cargo test` and `cargo clippy -- -D warnings`. Fix all failures.

## Notes

- From codebase review finding NEW-005.
- VariantType is in kam-call, not kam-core, so adding a variant is feasible
  without the "APPROVED" commit message requirement.
