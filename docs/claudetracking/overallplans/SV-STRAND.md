# SV-STRAND: SV-Aware Strand Bias Filter

**Status**: active
**Priority**: high
**Branch**: epic/SV-STRAND

## Goal

Prevent the StrandBias filter from removing correct inversion calls at moderate
and high VAF. At ≥2% VAF, the INV is reliably detected by the de Bruijn graph
but consistently removed by `StrandBias`. This makes inversions undetectable at
clinical VAF ranges even when signal is adequate.

## Background

The StrandBias filter applies Fisher's exact test on simplex molecule strand
counts: `(alt_fwd, alt_rev, ref_fwd, ref_rev)`. For SNVs and small indels, the
two alleles should be covered on both strands at a balanced ratio.

For inversion junction paths, this assumption fails by design:

- The de Bruijn graph walks the inversion path in one direction.
- Reads crossing the left junction are predominantly forward-oriented.
- Reads crossing the right junction are predominantly reverse-oriented.
- Very few reads span the entire 100bp inversion, so the graph accumulates
  most evidence from one strand end before the other.

The result: `min_simplex_fwd` or `min_simplex_rev` across the inversion path is
near 0 at the worst-covered k-mer, yielding a p-value near 0 and triggering
StrandBias for every INV call at ≥2% VAF.

Evidence (from `docs/research/sv_large_detection_investigation.md`):

| Truth VAF | NALT | Filter      |
|-----------|------|-------------|
| 2.0%      | 12   | StrandBias  |
| 5.0%      | 23   | StrandBias  |
| 10.0%     | 41   | StrandBias  |

All INV calls at ≥2% VAF are StrandBias-filtered regardless of NALT.

## Design

### Core change: SV-specific strand bias threshold

Add a `sv_strand_bias_threshold: f64` field to `CallerConfig` in
`kam-call/src/caller.rs`. Default: `1.0` (disabled — all INV/DUP/DEL calls
skip the strand bias filter).

The strand bias filter is not useful for SV paths because:
1. The path-level `min_simplex_fwd` / `min_simplex_rev` are not representative
   of strand balance at the variant site. They reflect the worst-covered k-mer,
   which is typically deep inside the inversion body.
2. The structural cause of strand imbalance (directional path walking) is not
   a sequencing artefact; it cannot be used to identify FPs.

Setting the default to 1.0 (i.e., all p-values pass) is appropriate because:
- SV false positives are already controlled by confidence and molecule count.
- No equivalent error pattern (a random 100bp inversion) exists in background
  sequencing noise.

### Caller change

In `assign_filter()`:

```rust
let eff_strand_bias_threshold = if is_sv {
    config.sv_strand_bias_threshold
} else {
    config.strand_bias_threshold
};
if strand_bias_p < eff_strand_bias_threshold {
    return VariantFilter::StrandBias;
}
```

### CLI change

Add `--sv-strand-bias-threshold <FLOAT>` to `CallerConfig` construction in
`kam/src/commands/call.rs` and `run.rs`. Default: `1.0` (disabled).

This is intentionally separate from `--strand-bias-threshold` so SNV/indel
behaviour is unchanged.

## Child tasks

| ID | File | Status |
|----|------|--------|
| SV-STRAND-001 | todo/sv_strand_001_config_field.md | todo |
| SV-STRAND-002 | todo/sv_strand_002_cli_wire.md | todo |
| SV-STRAND-003 | todo/sv_strand_003_benchmark_rerun.md | todo |

## Scope

- `kam-call/src/caller.rs` — add `sv_strand_bias_threshold` to `CallerConfig`,
  use in `assign_filter()`
- `kam/src/commands/call.rs` and `run.rs` — wire CLI flag
- `docs/benchmarking/sv/` — re-run and re-score SV suite

## Out of scope

- Changes to `kam-core` types
- Changes to path walking or k-mer scoring
- Junction-specific strand bias computation (future work)
