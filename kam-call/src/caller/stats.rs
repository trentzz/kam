//! Statistical helpers for variant calling: VAF estimation, strand bias, and posterior confidence.

use statrs::distribution::{Beta, ContinuousCDF};

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
    // Clamp k to m before subtraction. For SV types, k comes from
    // mean_variant_specific_molecules.round() while m = ref_k + k where
    // ref_k uses a different rounding path. Independent rounding can push k
    // above m, which would wrap on u32 subtraction.
    let k = k.min(m);
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
pub fn strand_bias_test(alt_fwd: u32, alt_rev: u32, ref_fwd: u32, ref_rev: u32) -> f64 {
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

/// Posterior probability that the variant is real (vs background error).
///
/// Uses a simple binomial likelihood ratio: probability of observing `k` or
/// more alt reads under signal (VAF = k/m) vs background error rate.
pub fn compute_confidence(k: u32, m: u32, background_error_rate: f64) -> f64 {
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

/// Natural log of the hypergeometric PMF for a 2×2 table.
///
/// `P(X = a) = C(row1,a) * C(row2, col1-a) / C(n, col1)`
pub(crate) fn hypergeometric_log_pmf(a: i64, row1: i64, col1: i64, n: i64) -> f64 {
    log_binom(row1, a) + log_binom(n - row1, col1 - a) - log_binom(n, col1)
}

/// Natural log of the binomial coefficient C(n, k).
fn log_binom(n: i64, k: i64) -> f64 {
    if k < 0 || k > n {
        return f64::NEG_INFINITY;
    }
    log_factorial(n) - log_factorial(k) - log_factorial(n - k)
}

/// Natural log of n! computed via the lgamma function (O(1)).
fn log_factorial(n: i64) -> f64 {
    if n <= 1 {
        return 0.0;
    }
    // ln(n!) = ln(Gamma(n+1)); statrs ln_gamma is O(1) via Lanczos approximation.
    statrs::function::gamma::ln_gamma(n as f64 + 1.0)
}

/// Log-likelihood of `k` successes in `m` Bernoulli trials with probability `p`.
fn log_binomial_likelihood(k: u32, m: u32, p: f64) -> f64 {
    if p <= 0.0 {
        if k == 0 {
            return 0.0;
        } else {
            return f64::NEG_INFINITY;
        }
    }
    if p >= 1.0 {
        if k == m {
            return 0.0;
        } else {
            return f64::NEG_INFINITY;
        }
    }
    let k_f = k as f64;
    let m_f = m as f64;
    k_f * p.ln() + (m_f - k_f) * (1.0 - p).ln()
}

// ─── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

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

    // Test: compute_confidence returns 0.0 when k=0.
    // Zero alt molecules provide zero evidence for the variant.
    #[test]
    fn compute_confidence_zero_alt_is_zero() {
        let c = compute_confidence(0, 1000, 1e-4);
        assert_eq!(c, 0.0, "zero alt molecules must give confidence 0.0");
    }

    // Test: compute_confidence returns 0.0 when m=0.
    // No data at all means no confidence.
    #[test]
    fn compute_confidence_zero_total_is_zero() {
        let c = compute_confidence(0, 0, 1e-4);
        assert_eq!(c, 0.0, "zero total molecules must give confidence 0.0");
    }

    // Test: strand_bias_test with all zeros returns 1.0.
    // When there are no reads at all, there can be no strand bias.
    #[test]
    fn strand_bias_all_zeros_returns_one() {
        let p = strand_bias_test(0, 0, 0, 0);
        assert_eq!(p, 1.0, "all-zero table must return p=1.0, got {p}");
    }

    // Test: estimate_vaf with zero total (m=0) returns (0, 0, 0).
    // This guards against division by zero when no molecules are observed.
    #[test]
    fn estimate_vaf_zero_total_returns_zeros() {
        let (vaf, lo, hi) = estimate_vaf(0, 0);
        assert_eq!(vaf, 0.0, "VAF must be 0 when total is 0");
        assert_eq!(lo, 0.0, "CI lower must be 0 when total is 0");
        assert_eq!(hi, 0.0, "CI upper must be 0 when total is 0");
    }

    // Test: estimate_vaf CI ordering invariant: lo < hi for any non-zero total.
    #[test]
    fn estimate_vaf_ci_ordering() {
        for (k, m) in [(1, 100), (50, 100), (99, 100), (0, 100), (5, 1000)] {
            let (_vaf, lo, hi) = estimate_vaf(k, m);
            assert!(
                lo < hi,
                "k={k}, m={m}: CI lower ({lo}) must be < CI upper ({hi})"
            );
            assert!(lo >= 0.0, "k={k}, m={m}: CI lower ({lo}) must be >= 0");
            assert!(hi <= 1.0, "k={k}, m={m}: CI upper ({hi}) must be <= 1.0");
        }
    }

    // Test (BUG-004): estimate_vaf must not panic or wrap when k > m.
    //
    // This can occur for SV types where k = mean_variant_specific_molecules.round()
    // and m = ref_k + k, but ref_k and k are rounded independently so k can
    // momentarily exceed m. The guard clamps k to m before subtraction.
    #[test]
    fn estimate_vaf_k_greater_than_m_clamps_not_panics() {
        // k=5 > m=3: without the guard this wraps on u32 subtraction.
        let (vaf, lo, hi) = estimate_vaf(5, 3);
        // Clamped: k becomes 3, vaf = 3/3 = 1.0.
        assert!((vaf - 1.0).abs() < 1e-9, "expected VAF=1.0, got {vaf}");
        // CI must still be valid floats in [0,1].
        assert!((0.0..=1.0).contains(&lo), "CI lower out of range: {lo}");
        assert!((0.0..=1.0).contains(&hi), "CI upper out of range: {hi}");
    }
}
