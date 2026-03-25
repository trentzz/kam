# CHEM-004: ReadStructure Presets

**Epic**: CHEM-CONFIG (docs/claudetracking/overallplans/CHEM-CONFIG.md)
**Priority**: critical
**Depends on**: chem_002_generalise_parser.md
**Status**: todo

## Goal

Add named constructor methods to `ReadStructure` (or `ChemistryConfig`) for
common sequencing chemistries, and wire the `ReadStructure` into the TOML
config chemistry section. Done looks like: `ReadStructure::twist_umi_duplex()`,
`ReadStructure::simplex_12bp()`, and `ReadStructure::simplex_9bp()` exist with
doc comments; the `chemistry.preset` TOML key (if set) initialises the correct
preset; and all three presets are tested.

## Steps

1. Read the `ReadStructure` or `ChemistryConfig` definition introduced in
   CHEM-002.

2. Add preset constructors:
   ```rust
   impl ReadStructure {
       /// Twist UMI duplex: 5 bp UMI, 2 bp skip (monotemplate spacer).
       pub fn twist_umi_duplex() -> Self { ... }

       /// Generic simplex protocol with 12 bp UMI and no skip.
       pub fn simplex_12bp() -> Self { ... }

       /// Generic simplex protocol with 9 bp UMI and no skip.
       pub fn simplex_9bp() -> Self { ... }
   }
   ```

3. Add an optional `preset` field to `ChemistryConfig` (in `kam/src/config.rs`
   if CONFIG-001 is done, or locally):
   ```toml
   [chemistry]
   preset = "twist-umi-duplex"  # or "simplex-12bp", "simplex-9bp"
   # umi_length and skip_length are ignored when preset is set
   ```

4. In `load_config` (or wherever `ChemistryConfig` is converted to
   `ReadStructure`), check the `preset` field first. If set, call the
   corresponding constructor. If not set, use `umi_length` and `skip_length`
   directly.

5. Add unit tests:
   - Each preset constructor produces the expected `umi_length` and
     `skip_length`.
   - TOML with `preset = "simplex-12bp"` produces a `ReadStructure` with
     `umi_length = 12`, `skip_length = 0`.
   - TOML with explicit `umi_length = 8` and no preset produces
     `ReadStructure { umi_length: 8, skip_length: 0 }`.

6. Add the `preset` field to both example config files created in CONFIG-004.

7. Run `cargo test` and `cargo clippy -- -D warnings`. Fix all failures.

## Notes

- If `preset` and explicit `umi_length` are both set, explicit fields win (or
  raise an error). Choose one behaviour and document it.
- The preset names use hyphens (`twist-umi-duplex`) in TOML for readability;
  the Rust constructors use underscores.
