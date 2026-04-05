# CHEM-001: Generalise UMI to Vec<u8> in kam-core (APPROVED)

**Epic**: CHEM-CONFIG (docs/claudetracking/overallplans/CHEM-CONFIG.md)
**Priority**: critical
**Depends on**: none
**Status**: todo

## Goal

Change `CanonicalUmiPair.umi_a`, `umi_b` and the equivalent fields in
`Molecule` from `[u8; 5]` to `Vec<u8>`. Update all code that constructs or
destructures these types. Done looks like: the codebase compiles, all tests
pass, and no `[u8; 5]` UMI type annotation remains in `kam-core`.

**NOTE**: This is an APPROVED core type change. The commit message must contain
"APPROVED".

## Steps

1. Read `kam-core/src/lib.rs` and all files in `kam-core/src/` to understand
   the full type surface.

2. Change the fields in `CanonicalUmiPair` and `Molecule`:
   ```rust
   // Before
   pub umi_a: [u8; 5],
   pub umi_b: [u8; 5],

   // After
   pub umi_a: Vec<u8>,
   pub umi_b: Vec<u8>,
   ```

3. Search for all construction sites of `CanonicalUmiPair` and `Molecule`
   across the workspace (use `cargo check` to find all type errors). Update
   each to use `Vec::from(...)` or collect from a slice.

4. Search for any code that pattern-matches on these fields as fixed-length
   arrays and update it to use slice operations.

5. Update any `PartialEq`, `Hash`, or serialisation impls that assumed fixed
   length. `Vec<u8>` derives these correctly, so explicit impls may be
   removable.

6. Run `cargo test` and `cargo clippy -- -D warnings`. Fix all failures.

## Notes

- `Vec<u8>` is heap-allocated and slightly slower than `[u8; 5]`. This is
  acceptable: UMI comparison is not on the hot path for k-mer indexing.
- Bincode serialisation of `Vec<u8>` produces a length-prefixed encoding,
  which differs from the fixed-array encoding. If serialised index files exist
  from prior runs, they will be incompatible. Document this in the commit
  message.
- The commit message must contain "APPROVED" because this modifies kam-core
  types.
