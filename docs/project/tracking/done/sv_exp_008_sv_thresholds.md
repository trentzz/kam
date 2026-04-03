# SV-EXP-008: Apply SV Thresholds to New Variant Types

**Epic**: SV-EXPAND (docs/claudetracking/overallplans/SV-EXPAND.md)
**Priority**: high
**Depends on**: sv_exp_001_new_variant_types.md
**Status**: todo

## Goal

Ensure `InvDel`, `NovelInsertion`, and `Fusion` are routed through the SV-specific
thresholds (`sv_min_confidence`, `sv_min_alt_molecules`, `sv_strand_bias_threshold`)
in `assign_filter`. Done looks like: the three new types use SV thresholds
in `assign_filter`; a unit test confirms each type is filtered by SV thresholds,
not SNV/indel thresholds.

## Steps

1. Read `kam-call/src/caller.rs`, focusing on `assign_filter` and the
   `is_sv_type` predicate.

2. Add `InvDel`, `NovelInsertion`, and `Fusion` to the `is_sv_type` match arm
   (or equivalent predicate). Verify `is_sv_type` now returns `true` for all
   five SV types (including existing `LargeDeletion`, `TandemDuplication`,
   `Inversion`).

3. Confirm that `assign_filter` uses `is_sv_type` to select the threshold set.
   No other changes should be needed.

4. Write unit tests:
   - `InvDel` with confidence 0.96 and `sv_min_confidence = 0.95` → PASS.
   - `InvDel` with confidence 0.96 and `min_confidence = 0.99` → would be
     filtered if SNV threshold applied, confirming SV threshold is used.
   - Same for `NovelInsertion` and `Fusion`.

5. Run `cargo test` and `cargo clippy -- -D warnings`. Fix all failures.

## Notes

- This task is intentionally small. `is_sv_type` is already implemented;
  this just adds three new variants to it.
- If `is_sv_type` is a method on `VariantType`, the addition in SV-EXP-001
  may have already added `todo!()` stubs. Replace those stubs here.
