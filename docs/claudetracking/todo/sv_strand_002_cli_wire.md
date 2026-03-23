# SV-STRAND-002: Wire sv_strand_bias_threshold CLI flag

**Epic**: SV-STRAND (docs/claudetracking/overallplans/SV-STRAND.md)
**Priority**: high
**Depends on**: SV-STRAND-001
**Status**: todo

## Goal

Expose `sv_strand_bias_threshold` as `--sv-strand-bias-threshold` in the
`call` and `run` subcommands. Done looks like: the flag is visible in `--help`,
defaults to 1.0, and is passed through to `CallerConfig` in both subcommands.

## Steps

1. Read `kam/src/cli.rs` and `kam/src/commands/call.rs` and `run.rs`.

2. In `cli.rs`, add to both `CallArgs` and `RunArgs`:
   ```rust
   /// Fisher p-value threshold for strand bias filter on SV-type variants.
   ///
   /// Defaults to 1.0 (disabled). Inversion junction reads are structurally
   /// strand-biased and the standard threshold is inappropriate for SV paths.
   /// Set to 0.0 to apply the same threshold as SNVs/indels.
   #[arg(long, default_value_t = 1.0f64)]
   pub sv_strand_bias_threshold: f64,
   ```

3. In `call.rs` and `run.rs`, where `CallerConfig` is constructed, add:
   ```rust
   sv_strand_bias_threshold: args.sv_strand_bias_threshold,
   ```

4. Update any test `RunArgs` or `CallArgs` struct initialisers to include
   `sv_strand_bias_threshold: 1.0`.

5. Run `cargo test` and `cargo clippy -- -D warnings`. Fix all failures.

## Notes

- Keep this separate from `--strand-bias-threshold` so SNV/indel behaviour
  is not accidentally changed.
- No benchmark re-run in this task. That is SV-STRAND-003.
