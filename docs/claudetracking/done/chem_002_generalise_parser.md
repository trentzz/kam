# CHEM-002: Generalise Parser to Configurable UMI/Skip Lengths

**Epic**: CHEM-CONFIG (docs/claudetracking/overallplans/CHEM-CONFIG.md)
**Priority**: critical
**Depends on**: chem_001_generalise_umi.md
**Status**: todo

## Goal

Change `ParsedReadPair` UMI and skip fields from fixed-length arrays to
`Vec<u8>`. Remove hardcoded length constants from the parser. The parser reads
`umi_length` and `skip_length` from a `ReadStructure` (or equivalent config
struct). Done looks like: the parser accepts a `ReadStructure` argument, all
length-specific code is removed, and a test with 12 bp UMI and 0 skip passes.

## Steps

1. Read `kam-assemble/src/parser.rs` in full.

2. Identify all hardcoded length constants (e.g. `const UMI_LEN: usize = 5`,
   `const SKIP_LEN: usize = 2`, or inline literals).

3. Add `ReadStructure` (or reuse `ChemistryConfig` from CONFIG-001) as a
   parameter to the parser. If CONFIG-001 is not yet done, define a minimal
   `ReadStructure { umi_length: usize, skip_length: usize }` locally and
   plan to replace it later.

4. Change `ParsedReadPair`:
   ```rust
   pub umi: Vec<u8>,
   pub skip: Vec<u8>,  // empty if skip_length == 0
   ```

5. Update parsing logic: slice bytes `[0..umi_length]` for UMI and
   `[umi_length..umi_length+skip_length]` for skip.

6. Remove any assertions of the form `assert_eq!(umi.len(), 5)`.

7. Update or add tests:
   - Existing Twist test (5 bp UMI, 2 bp skip) still passes.
   - New test: 12 bp UMI, 0 bp skip — correct UMI extraction.

8. Run `cargo test` and `cargo clippy -- -D warnings`. Fix all failures.

## Notes

- `ParsedReadPair` may not be a public type; check whether it is used outside
  `kam-assemble`. If internal only, the change is contained.
- The skip bytes are currently used for QC (monotemplate detection). Keep the
  skip field even if `skip_length = 0`; an empty `Vec<u8>` is fine.
