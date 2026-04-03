# BUG-004: u32 underflow in estimate_vaf if k > m

**Epic**: none (standalone bug fix)
**Priority**: medium
**Depends on**: none
**Status**: done

## Goal

Fix the potential u32 underflow in `estimate_vaf` in `kam-call/src/caller.rs`
where `(m - k) as f64` wraps if `k > m`. Done looks like: a guard clamps `k`
to `m` before the subtraction, with a test that exercises the edge case.

## Steps

1. Read `kam-call/src/caller.rs`, focusing on `estimate_vaf` (around line 327).

2. Add a guard before the computation:
   ```rust
   let k = k.min(m); // prevent underflow if rounding pushed k above m
   ```
   Add a brief comment explaining why this can occur (mean_variant_specific_molecules
   is rounded independently from the total).

3. Add a unit test: `estimate_vaf(k=5, m=3)` should not panic and should
   return a sensible result (VAF clamped to 1.0).

4. Move this task file to done and update Status.

5. Run `cargo test` and `cargo clippy -- -D warnings`. Fix all failures.

## Notes

- From codebase review finding NEW-004.
- This can occur with SV types where `k` comes from
  `mean_variant_specific_molecules.round()` and `m = ref_k + k` where
  `ref_k` uses a different metric.
