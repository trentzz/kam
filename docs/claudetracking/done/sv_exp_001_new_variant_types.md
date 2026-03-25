# SV-EXP-001: Add Fusion, InvDel, NovelInsertion to VariantType

**Epic**: SV-EXPAND (docs/claudetracking/overallplans/SV-EXPAND.md)
**Priority**: high
**Depends on**: none
**Status**: todo

## Goal

Add three new variants to the `VariantType` enum in `kam-core`: `Fusion`,
`InvDel`, and `NovelInsertion`. Update all match arms, `Display`, serialisation,
`is_sv` predicate, and filter logic. Done looks like: the three new variants
compile and round-trip through serialisation; `is_sv_type` returns `true` for
all three; no match is non-exhaustive.

## Steps

1. Read `kam-core/src/variant.rs` (or wherever `VariantType` is defined) in
   full.

2. Add the new variants:
   ```rust
   pub enum VariantType {
       // ... existing variants ...
       InvDel,
       NovelInsertion,
       Fusion,
   }
   ```

3. Update `Display` (or `fmt::Display` impl) to produce human-readable strings:
   - `InvDel` → `"InvDel"`
   - `NovelInsertion` → `"NovelInsertion"`
   - `Fusion` → `"Fusion"`

4. Update the `is_sv` (or `is_sv_type`) predicate to return `true` for all
   three new variants.

5. Search the workspace for all `match variant_type` blocks and add arms for
   the three new variants. For classification and output code that does not
   yet handle them, add `todo!()` stubs with a comment pointing to the
   implementing task.

6. Update VCF SVTYPE string mapping if it lives in `variant.rs`.

7. Run `cargo build` and `cargo clippy -- -D warnings`. Fix all failures.
   Do not run the full test suite yet (classification tests will fail until
   SV-EXP-002/003 are done).

## Notes

- Do not implement classification logic here. This task is structural only:
  add the variants, make the code compile, update all match arms.
- `serde` derive on `VariantType` will automatically handle the new variants.
  Check that serialised names match what downstream tools expect (use
  `#[serde(rename = "...")]` if needed).
