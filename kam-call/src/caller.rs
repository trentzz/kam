//! Statistical variant calling using binomial/Beta-posterior models.
//!
//! Given molecule-level evidence from reference and alternate paths, this
//! module estimates VAF with 95% credible intervals, tests for strand bias,
//! and applies configurable filters to produce calibrated variant calls.

use kam_pathfind::score::PathEvidence;
use statrs::distribution::{Beta, ContinuousCDF};

// ─── Public types ─────────────────────────────────────────────────────────────

/// A called variant with full statistical evidence.
///
/// # Example
/// ```
/// use kam_call::caller::{VariantFilter, VariantType};
/// let vt = VariantType::Snv;
/// let vf = VariantFilter::Pass;
/// assert_eq!(vt, VariantType::Snv);
/// assert_eq!(vf, VariantFilter::Pass);
/// ```
#[derive(Debug, Clone)]
pub struct VariantCall {
    /// Identifier for the target region.
    pub target_id: String,
    /// Variant class inferred from ref/alt sequence comparison.
    pub variant_type: VariantType,
    /// Reference allele sequence.
    pub ref_sequence: Vec<u8>,
    /// Alternate allele sequence.
    pub alt_sequence: Vec<u8>,
    /// VAF point estimate: `k / M`.
    pub vaf: f64,
    /// Lower bound of 95% Beta credible interval.
    pub vaf_ci_low: f64,
    /// Upper bound of 95% Beta credible interval.
    pub vaf_ci_high: f64,
    /// Number of molecules supporting the reference allele.
    pub n_molecules_ref: u32,
    /// Number of molecules supporting the alternate allele.
    pub n_molecules_alt: u32,
    /// Number of duplex molecules supporting the alternate allele.
    pub n_duplex_alt: u32,
    /// Number of simplex molecules supporting the alternate allele.
    pub n_simplex_alt: u32,
    /// Posterior probability that the variant is real.
    pub confidence: f64,
    /// Fisher's exact test p-value for strand bias.
    pub strand_bias_p: f64,
    /// Quality filter outcome.
    pub filter: VariantFilter,
}

/// Classification of a variant by allele length comparison.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VariantType {
    /// Single nucleotide variant.
    Snv,
    /// Insertion (alt longer than ref).
    Insertion,
    /// Deletion (alt shorter than ref).
    Deletion,
    /// Multi-nucleotide variant (same length, >1 difference).
    Mnv,
    /// Complex rearrangement.
    Complex,
}

/// Filter outcome for a variant call.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VariantFilter {
    /// Variant passes all filters.
    Pass,
    /// Significant strand bias detected.
    StrandBias,
    /// Posterior confidence below threshold.
    LowConfidence,
    /// Fewer duplex molecules than required.
    LowDuplex,
    /// UMI collision risk too high to trust.
    CollisionRisk,
}

/// Caller configuration with sensible defaults.
///
/// # Example
/// ```
/// use kam_call::caller::CallerConfig;
/// let cfg = CallerConfig::default();
/// assert_eq!(cfg.min_alt_molecules, 2);
/// ```
#[derive(Debug, Clone)]
pub struct CallerConfig {
    /// Minimum posterior probability to emit a `Pass` call. Default: 0.99.
    pub min_confidence: f64,
    /// Fisher p-value below which strand bias is flagged. Default: 0.01.
    pub strand_bias_threshold: f64,
    /// Minimum alt-supporting molecules required. Default: 2.
    pub min_alt_molecules: u32,
    /// Minimum duplex alt molecules required. Default: 0.
    pub min_alt_duplex: u32,
    /// Per-site background error rate used in the posterior. Default: 1e-4.
    pub background_error_rate: f64,
}

impl Default for CallerConfig {
    fn default() -> Self {
        Self {
            min_confidence: 0.99,
            strand_bias_threshold: 0.01,
            min_alt_molecules: 2,
            min_alt_duplex: 0,
            background_error_rate: 1e-4,
        }
    }
}

// ─── Public functions ─────────────────────────────────────────────────────────

/// Call a variant from reference and alternate path evidence.
///
/// Combines VAF estimation, strand-bias testing, posterior confidence, and
/// filter assignment into a single [`VariantCall`].
///
/// # Example
/// ```
/// use kam_call::caller::{call_variant, CallerConfig, VariantFilter};
/// use kam_pathfind::score::PathEvidence;
///
/// let ref_ev = PathEvidence {
///     min_molecules: 990, mean_molecules: 990.0,
///     min_duplex: 0, mean_duplex: 0.0,
///     total_simplex_fwd: 495, total_simplex_rev: 495,
///     mean_error_prob: 0.001,
/// };
/// let alt_ev = PathEvidence {
///     min_molecules: 10, mean_molecules: 10.0,
///     min_duplex: 5, mean_duplex: 5.0,
///     total_simplex_fwd: 5, total_simplex_rev: 5,
///     mean_error_prob: 0.001,
/// };
/// let call = call_variant("TP53", &ref_ev, &alt_ev, b"A", b"T", &CallerConfig::default());
/// assert_eq!(call.filter, VariantFilter::Pass);
/// ```
pub fn call_variant(
    target_id: &str,
    ref_evidence: &PathEvidence,
    alt_evidence: &PathEvidence,
    ref_seq: &[u8],
    alt_seq: &[u8],
    config: &CallerConfig,
) -> VariantCall {
    let k = alt_evidence.min_molecules;
    let total = ref_evidence.min_molecules + alt_evidence.min_molecules;

    let (vaf, vaf_ci_low, vaf_ci_high) = estimate_vaf(k, total);

    let strand_bias_p = strand_bias_test(
        alt_evidence.total_simplex_fwd,
        alt_evidence.total_simplex_rev,
        ref_evidence.total_simplex_fwd,
        ref_evidence.total_simplex_rev,
    );

    let confidence = compute_confidence(k, total, config.background_error_rate);

    let variant_type = classify_variant(ref_seq, alt_seq);

    let n_duplex_alt = alt_evidence.min_duplex;
    let n_simplex_alt = k.saturating_sub(n_duplex_alt);

    let filter = assign_filter(
        confidence,
        strand_bias_p,
        k,
        n_duplex_alt,
        config,
    );

    VariantCall {
        target_id: target_id.to_owned(),
        variant_type,
        ref_sequence: ref_seq.to_vec(),
        alt_sequence: alt_seq.to_vec(),
        vaf,
        vaf_ci_low,
        vaf_ci_high,
        n_molecules_ref: ref_evidence.min_molecules,
        n_molecules_alt: k,
        n_duplex_alt,
        n_simplex_alt,
        confidence,
        strand_bias_p,
        filter,
    }
}

/// Compute VAF point estimate and 95% Beta credible interval.
///
/// Uses the conjugate Beta posterior `Beta(k+1, M-k+1)` where `k` is the
/// number of alt-supporting molecules and `M` is the total molecule count.
/// Returns `(point_estimate, ci_lower_2.5%, ci_upper_97.5%)`.
///
/// # Example
/// ```
/// use kam_call::caller::estimate_vaf;
/// let (vaf, lo, hi) = estimate_vaf(5, 1000);
/// assert!((vaf - 0.005).abs() < 1e-6);
/// assert!(lo < 0.005 && hi > 0.005);
/// ```
pub fn estimate_vaf(k: u32, m: u32) -> (f64, f64, f64) {
    if m == 0 {
        return (0.0, 0.0, 0.0);
    }
    let point = k as f64 / m as f64;
    let alpha = k as f64 + 1.0;
    let beta_param = (m - k) as f64 + 1.0;
    // Beta distribution quantiles for the 95% credible interval.
    let dist = Beta::new(alpha, beta_param).expect("valid Beta parameters");
    let lo = dist.inverse_cdf(0.025);
    let hi = dist.inverse_cdf(0.975);
    (point, lo, hi)
}

/// Fisher's exact test for strand bias on a 2×2 contingency table.
///
/// Returns the two-tailed p-value computed by summing all hypergeometric
/// probability mass that is ≤ the probability of the observed table.
///
/// The table layout is:
///
/// ```text
/// alt_fwd  alt_rev
/// ref_fwd  ref_rev
/// ```
///
/// # Example
/// ```
/// use kam_call::caller::strand_bias_test;
/// // Perfectly balanced → high p-value (no bias).
/// let p = strand_bias_test(10, 10, 100, 100);
/// assert!(p > 0.5);
/// // All alt reads on one strand → extreme bias.
/// let p_biased = strand_bias_test(20, 0, 100, 100);
/// assert!(p_biased < 0.01);
/// ```
pub fn strand_bias_test(
    alt_fwd: u32,
    alt_rev: u32,
    ref_fwd: u32,
    ref_rev: u32,
) -> f64 {
    // 2×2 table:
    //   a = alt_fwd,  b = alt_rev   row1 = a+b
    //   c = ref_fwd,  d = ref_rev   row2 = c+d
    //   col1 = a+c,   col2 = b+d    n = a+b+c+d
    let a = alt_fwd as i64;
    let b = alt_rev as i64;
    let c = ref_fwd as i64;
    let d = ref_rev as i64;

    let row1 = a + b;
    let row2 = c + d;
    let col1 = a + c;
    let col2 = b + d;
    let n = row1 + row2;

    if n == 0 {
        return 1.0;
    }

    // Enumerate all possible values of cell `a` given the fixed margins.
    let a_min = 0_i64.max(row1 - col2);
    let a_max = row1.min(col1);

    // Log-probability of the observed table under the hypergeometric distribution.
    let log_p_obs = hypergeometric_log_pmf(a, row1, col1, n);

    // Sum probabilities of tables at least as extreme (p ≤ p_obs).
    let mut p_sum = 0.0_f64;
    let mut a_cur = a_min;
    while a_cur <= a_max {
        let log_p = hypergeometric_log_pmf(a_cur, row1, col1, n);
        if log_p <= log_p_obs + 1e-10 {
            p_sum += log_p.exp();
        }
        a_cur += 1;
    }

    p_sum.min(1.0)
}

/// Determine variant type by comparing reference and alternate sequences.
///
/// # Example
/// ```
/// use kam_call::caller::{classify_variant, VariantType};
/// assert_eq!(classify_variant(b"A", b"T"), VariantType::Snv);
/// assert_eq!(classify_variant(b"ACG", b"A"), VariantType::Deletion);
/// assert_eq!(classify_variant(b"A", b"ACG"), VariantType::Insertion);
/// ```
pub fn classify_variant(ref_seq: &[u8], alt_seq: &[u8]) -> VariantType {
    use std::cmp::Ordering;
    match ref_seq.len().cmp(&alt_seq.len()) {
        Ordering::Equal => {
            let diffs = ref_seq
                .iter()
                .zip(alt_seq.iter())
                .filter(|(r, a)| r != a)
                .count();
            match diffs {
                0 => VariantType::Snv, // identical — treat as SNV (no-op)
                1 => VariantType::Snv,
                _ => VariantType::Mnv,
            }
        }
        Ordering::Greater => VariantType::Deletion,
        Ordering::Less => VariantType::Insertion,
    }
}

// ─── Private helpers ──────────────────────────────────────────────────────────

/// Natural log of the hypergeometric PMF for a 2×2 table.
///
/// `P(X = a) = C(row1,a) * C(row2, col1-a) / C(n, col1)`
fn hypergeometric_log_pmf(a: i64, row1: i64, col1: i64, n: i64) -> f64 {
    log_binom(row1, a) + log_binom(n - row1, col1 - a) - log_binom(n, col1)
}

/// Natural log of the binomial coefficient C(n, k).
fn log_binom(n: i64, k: i64) -> f64 {
    if k < 0 || k > n {
        return f64::NEG_INFINITY;
    }
    log_factorial(n) - log_factorial(k) - log_factorial(n - k)
}

/// Natural log of n! computed via the sum of ln(i).
fn log_factorial(n: i64) -> f64 {
    if n <= 1 {
        return 0.0;
    }
    (2..=n).map(|i| (i as f64).ln()).sum()
}

/// Posterior probability that the variant is real (vs background error).
///
/// Uses a simple binomial likelihood ratio: probability of observing `k` or
/// more alt reads under signal (VAF = k/m) vs background error rate.
fn compute_confidence(k: u32, m: u32, background_error_rate: f64) -> f64 {
    if k == 0 || m == 0 {
        return 0.0;
    }

    let vaf = k as f64 / m as f64;

    // Compute log-likelihood of observing k successes out of m under VAF signal.
    let log_l_signal = log_binomial_likelihood(k, m, vaf);

    // Compute log-likelihood under background error rate.
    let log_l_background = log_binomial_likelihood(k, m, background_error_rate);

    // Bayes factor (equal priors → posterior ∝ likelihood ratio).
    let log_ratio = log_l_signal - log_l_background;

    // Convert to posterior probability with equal priors.
    let exp_ratio = log_ratio.exp();
    if exp_ratio.is_infinite() {
        return 1.0;
    }
    exp_ratio / (1.0 + exp_ratio)
}

/// Log-likelihood of `k` successes in `m` Bernoulli trials with probability `p`.
fn log_binomial_likelihood(k: u32, m: u32, p: f64) -> f64 {
    if p <= 0.0 {
        if k == 0 { return 0.0; } else { return f64::NEG_INFINITY; }
    }
    if p >= 1.0 {
        if k == m { return 0.0; } else { return f64::NEG_INFINITY; }
    }
    let k_f = k as f64;
    let m_f = m as f64;
    k_f * p.ln() + (m_f - k_f) * (1.0 - p).ln()
}

/// Assign a filter label given statistical evidence and thresholds.
fn assign_filter(
    confidence: f64,
    strand_bias_p: f64,
    n_alt: u32,
    n_duplex_alt: u32,
    config: &CallerConfig,
) -> VariantFilter {
    if n_alt < config.min_alt_molecules {
        return VariantFilter::LowConfidence;
    }
    if n_duplex_alt < config.min_alt_duplex {
        return VariantFilter::LowDuplex;
    }
    if strand_bias_p < config.strand_bias_threshold {
        return VariantFilter::StrandBias;
    }
    if confidence < config.min_confidence {
        return VariantFilter::LowConfidence;
    }
    VariantFilter::Pass
}

// ─── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use kam_pathfind::score::PathEvidence;

    fn make_path_evidence(
        min_molecules: u32,
        min_duplex: u32,
        fwd: u32,
        rev: u32,
    ) -> PathEvidence {
        PathEvidence {
            min_molecules,
            mean_molecules: min_molecules as f32,
            min_duplex,
            mean_duplex: min_duplex as f32,
            total_simplex_fwd: fwd,
            total_simplex_rev: rev,
            mean_error_prob: 0.001,
        }
    }

    // Test 1: estimate_vaf(5, 1000) → VAF ≈ 0.005, CI contains 0.005.
    #[test]
    fn estimate_vaf_small_fraction() {
        let (vaf, lo, hi) = estimate_vaf(5, 1000);
        assert!((vaf - 0.005).abs() < 1e-6, "VAF should be 0.005, got {vaf}");
        assert!(lo < 0.005, "CI lower {lo} should be below 0.005");
        assert!(hi > 0.005, "CI upper {hi} should be above 0.005");
    }

    // Test 2: estimate_vaf(0, 1000) → VAF = 0, CI lower near 0.
    //
    // With Beta(1, 1001) the 2.5th percentile is very small (~2.5e-5) but not
    // exactly zero.  The spec says "CI lower = 0" meaning it should be close
    // to zero rather than strictly 0.0.
    #[test]
    fn estimate_vaf_zero_alt() {
        let (vaf, lo, _hi) = estimate_vaf(0, 1000);
        assert_eq!(vaf, 0.0, "VAF should be 0");
        assert!(lo < 1e-3, "CI lower should be near 0, got {lo}");
    }

    // Test 3: strand_bias_test with balanced strands → high p-value.
    #[test]
    fn strand_bias_balanced_is_nonsignificant() {
        let p = strand_bias_test(10, 10, 100, 100);
        assert!(p > 0.5, "balanced strands should yield p > 0.5, got {p}");
    }

    // Test 4: strand_bias_test with all-one-strand → low p-value.
    #[test]
    fn strand_bias_all_one_strand_is_significant() {
        let p = strand_bias_test(20, 0, 100, 100);
        assert!(p < 0.01, "all-one-strand should yield p < 0.01, got {p}");
    }

    // Test 5: classify_variant with same-length, 1 diff → SNV.
    #[test]
    fn classify_snv() {
        assert_eq!(classify_variant(b"A", b"T"), VariantType::Snv);
        assert_eq!(classify_variant(b"ACGT", b"ACTT"), VariantType::Snv);
    }

    // Test 6: classify_variant alt shorter → Deletion.
    #[test]
    fn classify_deletion() {
        assert_eq!(classify_variant(b"ACG", b"A"), VariantType::Deletion);
    }

    // Test 7: classify_variant alt longer → Insertion.
    #[test]
    fn classify_insertion() {
        assert_eq!(classify_variant(b"A", b"ACG"), VariantType::Insertion);
    }

    // Test 8: call_variant with strong evidence → PASS filter.
    #[test]
    fn call_variant_strong_evidence_passes() {
        let ref_ev = make_path_evidence(990, 0, 495, 495);
        let alt_ev = make_path_evidence(10, 5, 5, 5);
        let cfg = CallerConfig::default();
        let call = call_variant("TP53", &ref_ev, &alt_ev, b"A", b"T", &cfg);
        assert_eq!(call.filter, VariantFilter::Pass);
        assert!((call.vaf - 10.0 / 1000.0).abs() < 1e-6);
    }

    // Test 9: call_variant with strand bias → StrandBias filter.
    #[test]
    fn call_variant_strand_bias_flagged() {
        // Alt reads are all on the forward strand.
        let ref_ev = make_path_evidence(990, 0, 495, 495);
        let alt_ev = PathEvidence {
            min_molecules: 10,
            mean_molecules: 10.0,
            min_duplex: 5,
            mean_duplex: 5.0,
            total_simplex_fwd: 10, // all forward
            total_simplex_rev: 0,
            mean_error_prob: 0.001,
        };
        let cfg = CallerConfig::default();
        let call = call_variant("TP53", &ref_ev, &alt_ev, b"A", b"T", &cfg);
        assert_eq!(call.filter, VariantFilter::StrandBias);
    }

    // Test 10: call_variant with low molecule count → LowConfidence filter.
    #[test]
    fn call_variant_low_molecules_flagged() {
        // Only 1 alt molecule — below the default minimum of 2.
        let ref_ev = make_path_evidence(999, 0, 500, 499);
        let alt_ev = make_path_evidence(1, 0, 1, 0);
        let cfg = CallerConfig::default();
        let call = call_variant("TP53", &ref_ev, &alt_ev, b"A", b"T", &cfg);
        assert_eq!(call.filter, VariantFilter::LowConfidence);
    }
}
