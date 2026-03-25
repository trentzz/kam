# CONFIG-004: Example Config Files

**Epic**: CONFIG-TOML (docs/claudetracking/overallplans/CONFIG-TOML.md)
**Priority**: critical
**Depends on**: config_001_toml_schema.md
**Status**: todo

## Goal

Create two fully documented example config files in `examples/`. Done looks
like: `examples/twist-umi-duplex.toml` and `examples/simplex-umi-12bp.toml`
exist, every field has an inline comment explaining what it does and its valid
range, and `kam run --config examples/twist-umi-duplex.toml --r1 ... --r2 ...`
runs without error.

## Steps

1. Create `examples/twist-umi-duplex.toml`:
   - All sections present: `[chemistry]`, `[assembly]`, `[indexing]`,
     `[calling]`, `[output]`, `[input]`.
   - Values set to the Twist UMI duplex defaults.
   - Every field commented with: purpose, valid range, and consequence of
     changing it.
   - Example:
     ```toml
     [chemistry]
     # Number of UMI bases at the start of each read (Twist = 5).
     umi_length = 5
     # Number of spacer bases after the UMI (Twist monotemplate spacer = 2).
     skip_length = 2
     ```

2. Create `examples/simplex-umi-12bp.toml`:
   - Same structure, values for a simplex 12 bp UMI protocol with no skip.
   - Note which fields differ from Twist defaults and why.

3. Verify both files parse without error:
   ```bash
   toml-check examples/twist-umi-duplex.toml  # or use a small Rust test
   ```
   Alternatively, write a `#[test]` in `config.rs` that deserialises each file.

4. Add a reference to both files in `docs/manual/` (or `README.md` if no
   manual exists yet).

## Notes

- The `[input]` section in example files should leave paths empty or use
  placeholder strings like `""` so users know to fill them in.
- Comments in TOML use `#`. Use a blank line between sections for readability.
- This task can be done in parallel with CONFIG-002 and CONFIG-003 since it
  only depends on the schema (CONFIG-001).
