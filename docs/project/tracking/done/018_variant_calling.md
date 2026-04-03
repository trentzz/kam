# Task 018: Statistical Variant Calling

## Location
`kam-call/src/caller.rs`

## What to implement

Binomial model for variant calling: given scored paths (reference and variant), estimate VAF with credible intervals, test against background, compute confidence.

See `docs/research/statistical_calling_models.md` for the math.

## Interface

```rust
/// A called variant with statistical evidence
#[derive(Debug, Clone)]
pub struct VariantCall {
    pub target_id: String,
    pub variant_type: VariantType,
    pub ref_sequence: Vec<u8>,
    pub alt_sequence: Vec<u8>,
    pub vaf: f64,                    // point estimate (k / M)
    pub vaf_ci_low: f64,            // 95% credible interval lower bound
    pub vaf_ci_high: f64,           // 95% credible interval upper bound
    pub n_molecules_ref: u32,
    pub n_molecules_alt: u32,
    pub n_duplex_alt: u32,
    pub n_simplex_alt: u32,
    pub confidence: f64,             // posterior P(variant is real)
    pub strand_bias_p: f64,          // Fisher's exact test p-value
    pub filter: VariantFilter,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VariantType {
    Snv,
    Insertion,
    Deletion,
    Mnv,
    Complex,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VariantFilter {
    Pass,
    StrandBias,
    LowConfidence,
    LowDuplex,
    CollisionRisk,
}

/// Configuration for the caller
#[derive(Debug, Clone)]
pub struct CallerConfig {
    pub min_confidence: f64,         // minimum posterior prob to call PASS (default: 0.99)
    pub strand_bias_threshold: f64,  // p-value below which to flag (default: 0.01)
    pub min_alt_molecules: u32,      // minimum molecules supporting alt (default: 2)
    pub min_alt_duplex: u32,         // minimum duplex molecules for alt (default: 0)
    pub background_error_rate: f64,  // per-site background (default: 1e-4)
}

impl Default for CallerConfig { ... }

/// Call a variant from reference and alt path evidence.
pub fn call_variant(
    target_id: &str,
    ref_evidence: &PathEvidence,
    alt_evidence: &PathEvidence,
    ref_seq: &[u8],
    alt_seq: &[u8],
    config: &CallerConfig,
) -> VariantCall;

/// Compute VAF point estimate and 95% credible interval using Beta posterior.
/// k = alt molecules, M = total molecules
pub fn estimate_vaf(k: u32, m: u32) -> (f64, f64, f64);

/// Fisher's exact test for strand bias (2x2 table).
/// Returns p-value.
pub fn strand_bias_test(
    alt_fwd: u32, alt_rev: u32,
    ref_fwd: u32, ref_rev: u32,
) -> f64;

/// Determine variant type by comparing ref and alt sequences
pub fn classify_variant(ref_seq: &[u8], alt_seq: &[u8]) -> VariantType;
```

## Tests required

1. estimate_vaf(5, 1000) → VAF ≈ 0.005, CI contains 0.005
2. estimate_vaf(0, 1000) → VAF = 0, CI lower = 0
3. strand_bias_test with balanced strands → high p-value (no bias)
4. strand_bias_test with all-one-strand → low p-value
5. classify_variant: same length, 1 diff → SNV
6. classify_variant: alt shorter → Deletion
7. classify_variant: alt longer → Insertion
8. call_variant with strong evidence → PASS filter
9. call_variant with strand bias → StrandBias filter
10. call_variant with low molecule count → LowConfidence filter

## Notes

For Fisher's exact test, implement a simple version using the hypergeometric distribution. The `statrs` crate has `Hypergeometric` distribution. Alternatively, for a 2x2 table, compute the exact p-value by summing probabilities of tables at least as extreme.

For the Beta posterior: `Beta(k+1, M-k+1)`, use `statrs` crate's `Beta` distribution for quantiles.

## Definition of done

- `cargo test -p kam-call` passes
- `cargo clippy -p kam-call -- -D warnings` passes
- Doc comments with examples
- Add `statrs` to kam-call dependencies
