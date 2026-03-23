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
    /// Insertion (alt longer than ref, < 50 bp).
    Insertion,
    /// Deletion (alt shorter than ref, < 50 bp).
    Deletion,
    /// Multi-nucleotide variant (same length, >1 difference).
    Mnv,
    /// Complex rearrangement.
    Complex,
    /// Large deletion (alt shorter than ref by ≥ 50 bp).
    LargeDeletion,
    /// Tandem duplication (alt longer than ref by ≥ 50 bp).
    TandemDuplication,
    /// Inversion (same length, alt is reverse-complement of ref segment).
    Inversion,
}

/// Minimum indel length to classify as a structural variant.
pub const SV_LENGTH_THRESHOLD: usize = 50;

/// Filter outcome for a variant call.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VariantFilter {
    /// Variant passes all filters.
    Pass,
    /// Significant strand bias detected.
    StrandBias,
    /// Posterior confidence below threshold.
    LowConfidence,
    /// Fewer duplex molecules than required at the variant site.
    LowDuplex,
    /// UMI collision risk too high to trust.
    CollisionRisk,
    /// VAF exceeds the configured maximum (likely germline).
    HighVaf,
    /// Variant does not match any entry in the target variants set.
    ///
    /// Applied in tumour-informed monitoring mode (`--target-variants`).
    /// A call at this locus passed all statistical filters but the
    /// called allele is not one of the expected somatic variants.
    NotTargeted,
}

/// Caller configuration with sensible defaults.
///
/// # Example
/// ```
/// use kam_call::caller::CallerConfig;
/// let cfg = CallerConfig::default();
/// assert_eq!(cfg.min_alt_molecules, 2);
/// assert_eq!(cfg.min_alt_duplex, 0);
/// assert_eq!(cfg.max_vaf, None);
/// assert_eq!(cfg.sv_min_confidence, 0.95);
/// assert_eq!(cfg.sv_min_alt_molecules, 1);
/// ```
#[derive(Debug, Clone)]
pub struct CallerConfig {
    /// Minimum posterior probability to emit a `Pass` call. Default: 0.99.
    pub min_confidence: f64,
    /// Fisher p-value below which strand bias is flagged. Default: 0.01.
    pub strand_bias_threshold: f64,
    /// Minimum alt-supporting molecules required. Default: 2.
    ///
    /// If `min_alt_duplex_for_single` is set (not None) and the variant has at
    /// least that many duplex molecules at the variant site, a single-molecule
    /// call (n_alt = 1) is accepted.  Otherwise this threshold applies strictly.
    pub min_alt_molecules: u32,
    /// Minimum variant-specific duplex molecules required.
    ///
    /// Uses the minimum duplex count at variant-specific k-mers only (k-mers
    /// in the alt path that are absent from the reference path).  This is a
    /// calibrated signal: the old `mean_duplex` field included anchor k-mers
    /// covered by many reference molecules and was always ≥ 1, making the
    /// filter ineffective.  Default: 1.
    pub min_alt_duplex: u32,
    /// Minimum variant-specific duplex to accept a single-molecule call.
    ///
    /// When set, a call with n_alt = 1 passes the molecule threshold if its
    /// variant-specific duplex count is at least this value.  Duplex
    /// confirmation on a single molecule is strong evidence that the variant
    /// is real rather than a sequencing error on one strand.
    /// `None` disables single-molecule calls regardless of duplex support.
    /// Default: `Some(1)`.
    pub min_alt_duplex_for_single: Option<u32>,
    /// Maximum VAF for a `Pass` call.  Calls above this are filtered as
    /// `HighVaf`.  Useful for somatic variant calling in ctDNA where germline
    /// heterozygous variants (VAF ≈ 0.5) should be excluded.
    /// `None` disables the filter.  Default: `None`.
    pub max_vaf: Option<f64>,
    /// Per-site background error rate used in the posterior. Default: 1e-4.
    pub background_error_rate: f64,
    /// Minimum posterior probability for structural variant types
    /// (`LargeDeletion`, `TandemDuplication`, `Inversion`). Default: 0.95.
    ///
    /// A large structural event with 2 supporting molecules is qualitatively
    /// stronger evidence than 2 molecules supporting a single-base change —
    /// the probability that a random sequencing error mimics a 100 bp deletion
    /// or inversion is negligible. A lower threshold than `min_confidence` is
    /// appropriate for SV types.
    pub sv_min_confidence: f64,
    /// Minimum alt molecule count for structural variant types. Default: 1.
    ///
    /// SVs with even a single supporting molecule can be reported in
    /// monitoring mode where the target allele is pre-specified. In discovery
    /// mode, the confidence filter still governs whether the call is PASS.
    pub sv_min_alt_molecules: u32,
    /// Fisher p-value threshold for strand bias filter on SV-type variants.
    ///
    /// Defaults to 1.0 (disabled). Inversion junction reads are structurally
    /// strand-biased due to directional path walking in the de Bruijn graph.
    /// The standard `strand_bias_threshold` is inappropriate for SV paths.
    ///
    /// Set to 0.0 to apply the same threshold as SNVs/indels.
    pub sv_strand_bias_threshold: f64,
}

impl Default for CallerConfig {
    fn default() -> Self {
        Self {
            min_confidence: 0.99,
            strand_bias_threshold: 0.01,
            min_alt_molecules: 2,
            // 0 disables the LowDuplex hard filter.  At typical depths (2M reads,
            // 6.6% duplex rate) roughly 25–36% of genuine 2% VAF variants have
            // zero duplex coverage at the variant site, so a threshold of 1 would
            // remove too many true positives.  Enable by setting to 1 when running
            // at higher duplex rates (≥15%) or deeper sequencing (≥5M reads).
            min_alt_duplex: 0,
            min_alt_duplex_for_single: Some(1),
            max_vaf: None,
            background_error_rate: 1e-4,
            sv_min_confidence: 0.95,
            sv_min_alt_molecules: 1,
            sv_strand_bias_threshold: 1.0,
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
///     min_variant_specific_duplex: 0,
///     mean_variant_specific_molecules: 0.0,
///     min_simplex_fwd: 495, min_simplex_rev: 495,
///     mean_error_prob: 0.001,
/// };
/// let alt_ev = PathEvidence {
///     min_molecules: 10, mean_molecules: 10.0,
///     min_duplex: 5, mean_duplex: 5.0,
///     min_variant_specific_duplex: 5,
///     mean_variant_specific_molecules: 10.0,
///     min_simplex_fwd: 5, min_simplex_rev: 5,
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
    // Identical ref and alt: no variant to call. Return a filtered no-op call
    // immediately so downstream code never sees a spurious Snv at 100% VAF.
    if ref_seq == alt_seq {
        return VariantCall {
            target_id: target_id.to_owned(),
            variant_type: VariantType::Snv,
            ref_sequence: ref_seq.to_vec(),
            alt_sequence: alt_seq.to_vec(),
            vaf: 0.0,
            vaf_ci_low: 0.0,
            vaf_ci_high: 0.0,
            n_molecules_ref: ref_evidence.min_molecules,
            n_molecules_alt: 0,
            n_duplex_alt: 0,
            n_simplex_alt: 0,
            confidence: 0.0,
            strand_bias_p: 1.0,
            filter: VariantFilter::LowConfidence,
        };
    }

    // Classify the variant first so we can choose the appropriate evidence metric.
    let variant_type = classify_variant(ref_seq, alt_seq);

    // SV types use mean_variant_specific_molecules as the alt molecule count.
    // min_molecules bottlenecks at 1–3 for long SV paths (70–150 k-mers) even
    // at moderate VAF because only a handful of variant-specific k-mers are
    // observed at each molecule.  The mean over variant-specific k-mers gives
    // a calibrated count that is robust to the length of the SV path.
    let is_sv = matches!(
        variant_type,
        VariantType::LargeDeletion | VariantType::TandemDuplication | VariantType::Inversion
    );
    let k = if is_sv {
        alt_evidence.mean_variant_specific_molecules.round() as u32
    } else {
        alt_evidence.min_molecules
    };
    let total = if is_sv {
        // Use mean_molecules for the reference: the reference path spans the
        // whole target window and does not have the min bottleneck.
        let ref_k = ref_evidence.mean_molecules.round() as u32;
        ref_k + k
    } else {
        ref_evidence.min_molecules + alt_evidence.min_molecules
    };

    let (vaf, vaf_ci_low, vaf_ci_high) = estimate_vaf(k, total);

    let strand_bias_p = strand_bias_test(
        alt_evidence.min_simplex_fwd,
        alt_evidence.min_simplex_rev,
        ref_evidence.min_simplex_fwd,
        ref_evidence.min_simplex_rev,
    );

    let confidence = compute_confidence(k, total, config.background_error_rate);

    // Use the minimum duplex count at variant-specific k-mers only.
    // The old mean_duplex across all ~70 path k-mers was dominated by anchor
    // k-mers covered by reference molecules, making it always ≥ 1 and
    // rendering the LowDuplex filter ineffective.  min_variant_specific_duplex
    // reflects actual duplex support at the variant site.
    let n_duplex_alt = alt_evidence.min_variant_specific_duplex;
    let n_simplex_alt = k.saturating_sub(n_duplex_alt);

    let filter = assign_filter(
        confidence,
        strand_bias_p,
        k,
        n_duplex_alt,
        vaf,
        variant_type,
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

/// Determine variant type by comparing reference and alternate sequences.
///
/// Deletions and insertions of 50 bp or more are classified as structural
/// variants (`LargeDeletion` and `TandemDuplication` respectively). Same-length
/// variants where the alt is the reverse complement of the ref are classified as
/// `Inversion`.
///
/// # Example
/// ```
/// use kam_call::caller::{classify_variant, VariantType};
/// assert_eq!(classify_variant(b"A", b"T"), VariantType::Snv);
/// assert_eq!(classify_variant(b"ACG", b"A"), VariantType::Deletion);
/// assert_eq!(classify_variant(b"A", b"ACG"), VariantType::Insertion);
/// // Large structural variants:
/// let long_ref: Vec<u8> = b"ACGT".repeat(25).to_vec(); // 100 bp
/// let short_alt: Vec<u8> = b"ACGT".repeat(10).to_vec(); // 40 bp — 60 bp deletion
/// assert_eq!(classify_variant(&long_ref, &short_alt), VariantType::LargeDeletion);
/// ```
pub fn classify_variant(ref_seq: &[u8], alt_seq: &[u8]) -> VariantType {
    use std::cmp::Ordering;
    match ref_seq.len().cmp(&alt_seq.len()) {
        Ordering::Equal => {
            // Check for inversion: alt is reverse complement of ref.
            // Full-path case: entire window is RC (rare).
            if is_reverse_complement(ref_seq, alt_seq) {
                return VariantType::Inversion;
            }
            // Partial-inversion case: a contiguous central segment is RC'd
            // while the flanking bases match. This is the common case for
            // targeted SV detection where the target window is wider than
            // the inverted region.
            if let Some(inv_len) = partial_inversion_len(ref_seq, alt_seq) {
                if inv_len >= SV_LENGTH_THRESHOLD {
                    return VariantType::Inversion;
                }
            }
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
        Ordering::Greater => {
            if ref_seq.len() - alt_seq.len() >= SV_LENGTH_THRESHOLD {
                VariantType::LargeDeletion
            } else {
                VariantType::Deletion
            }
        }
        Ordering::Less => {
            if alt_seq.len() - ref_seq.len() >= SV_LENGTH_THRESHOLD {
                VariantType::TandemDuplication
            } else {
                VariantType::Insertion
            }
        }
    }
}

/// Detect a partial inversion: a contiguous central segment where alt is the
/// reverse complement of the corresponding ref region, with matching flanks.
///
/// Returns `Some(length)` where `length` is the number of bases in the inverted
/// segment if such a segment exists, or `None` otherwise.
///
/// Only considers segments of length ≥ 2 (required by `is_reverse_complement`).
fn partial_inversion_len(ref_seq: &[u8], alt_seq: &[u8]) -> Option<usize> {
    if ref_seq.len() != alt_seq.len() {
        return None;
    }
    let left = ref_seq
        .iter()
        .zip(alt_seq.iter())
        .position(|(r, a)| r != a)?;
    let right = ref_seq
        .iter()
        .zip(alt_seq.iter())
        .rposition(|(r, a)| r != a)?;
    if right < left {
        return None;
    }
    if is_reverse_complement(&ref_seq[left..=right], &alt_seq[left..=right]) {
        Some(right - left + 1)
    } else {
        None
    }
}

/// Return true if `alt` is the reverse complement of `ref`.
///
/// Only returns true for sequences of length ≥ 2 where the entire alt is the
/// reverse complement of the entire ref. This detects inversion paths in the de
/// Bruijn graph where both strands of a segment are traversed.
fn is_reverse_complement(ref_seq: &[u8], alt_seq: &[u8]) -> bool {
    if ref_seq.len() < 2 || ref_seq.len() != alt_seq.len() {
        return false;
    }
    ref_seq
        .iter()
        .zip(alt_seq.iter().rev())
        .all(|(&r, &a)| complement(r) == a)
}

/// DNA complement (ACGT only; other bytes return themselves unchanged).
#[inline]
fn complement(b: u8) -> u8 {
    match b {
        b'A' | b'a' => b'T',
        b'T' | b't' => b'A',
        b'C' | b'c' => b'G',
        b'G' | b'g' => b'C',
        other => other,
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

/// Natural log of n! computed via the lgamma function (O(1)).
fn log_factorial(n: i64) -> f64 {
    if n <= 1 {
        return 0.0;
    }
    // ln(n!) = ln(Gamma(n+1)); statrs ln_gamma is O(1) via Lanczos approximation.
    statrs::function::gamma::ln_gamma(n as f64 + 1.0)
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

/// Assign a filter label given statistical evidence and thresholds.
///
/// SV types (`LargeDeletion`, `TandemDuplication`, `Inversion`) use
/// `sv_min_confidence` and `sv_min_alt_molecules` instead of the SNV/indel
/// defaults. A large structural event requires fewer supporting molecules to
/// reach the same confidence as a single-base change, because background
/// sequencing errors cannot produce a 50–200 bp structural signature.
fn assign_filter(
    confidence: f64,
    strand_bias_p: f64,
    n_alt: u32,
    n_duplex_alt: u32,
    vaf: f64,
    variant_type: VariantType,
    config: &CallerConfig,
) -> VariantFilter {
    let is_sv = matches!(
        variant_type,
        VariantType::LargeDeletion | VariantType::TandemDuplication | VariantType::Inversion
    );
    let eff_min_molecules = if is_sv {
        config.sv_min_alt_molecules
    } else {
        config.min_alt_molecules
    };
    let eff_min_confidence = if is_sv {
        config.sv_min_confidence
    } else {
        config.min_confidence
    };

    // Molecule count: accept single-molecule calls when duplex-confirmed.
    let min_mols_ok = n_alt >= eff_min_molecules
        || (n_alt == 1
            && config
                .min_alt_duplex_for_single
                .is_some_and(|min_d| n_duplex_alt >= min_d));
    if !min_mols_ok {
        return VariantFilter::LowConfidence;
    }
    if n_duplex_alt < config.min_alt_duplex {
        return VariantFilter::LowDuplex;
    }
    let eff_strand_bias_threshold = if is_sv {
        config.sv_strand_bias_threshold
    } else {
        config.strand_bias_threshold
    };
    // A threshold of 1.0 disables the strand-bias filter (all p-values are < 1.0).
    if eff_strand_bias_threshold < 1.0 && strand_bias_p < eff_strand_bias_threshold {
        return VariantFilter::StrandBias;
    }
    if confidence < eff_min_confidence {
        return VariantFilter::LowConfidence;
    }
    if let Some(max_vaf) = config.max_vaf {
        if vaf > max_vaf {
            return VariantFilter::HighVaf;
        }
    }
    VariantFilter::Pass
}

// ─── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use kam_pathfind::score::PathEvidence;

    fn make_path_evidence(min_molecules: u32, min_duplex: u32, fwd: u32, rev: u32) -> PathEvidence {
        PathEvidence {
            min_molecules,
            mean_molecules: min_molecules as f32,
            min_duplex,
            mean_duplex: min_duplex as f32,
            // For tests, set variant-specific duplex equal to min_duplex so
            // that the filter behaves the same as the old mean_duplex logic.
            min_variant_specific_duplex: min_duplex,
            // For tests, set variant-specific mean equal to the overall mean.
            mean_variant_specific_molecules: min_molecules as f32,
            min_simplex_fwd: fwd,
            min_simplex_rev: rev,
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
            min_variant_specific_duplex: 5,
            mean_variant_specific_molecules: 10.0,
            min_simplex_fwd: 10, // all forward
            min_simplex_rev: 0,
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
        // The min_alt_molecules check fires before min_alt_duplex.
        let ref_ev = make_path_evidence(999, 0, 500, 499);
        let alt_ev = make_path_evidence(1, 0, 1, 0);
        let cfg = CallerConfig::default();
        let call = call_variant("TP53", &ref_ev, &alt_ev, b"A", b"T", &cfg);
        assert_eq!(call.filter, VariantFilter::LowConfidence);
    }

    // Test 11: call_variant with min_alt_duplex=1 and no variant-specific duplex → LowDuplex.
    //
    // The default is min_alt_duplex=0 (disabled); this test uses an explicit
    // config with threshold=1 to verify the filter works when enabled.
    #[test]
    fn call_variant_no_duplex_flagged_when_threshold_set() {
        let ref_ev = make_path_evidence(990, 0, 495, 495);
        let alt_ev = make_path_evidence(10, 0, 5, 5); // min_variant_specific_duplex = 0
        let cfg = CallerConfig {
            min_alt_duplex: 1,
            ..CallerConfig::default()
        };
        let call = call_variant("TP53", &ref_ev, &alt_ev, b"A", b"T", &cfg);
        assert_eq!(call.filter, VariantFilter::LowDuplex);
    }

    // Test 12: single molecule with duplex confirmation passes when min_alt_duplex_for_single is set.
    //
    // Uses a low-coverage target (13 ref + 1 alt = 14 molecules, VAF ≈ 7%),
    // analogous to the EGFR T790M case in the FN investigation.  At this depth
    // the posterior confidence for n_alt=1 exceeds 0.99, so the call passes.
    #[test]
    fn single_molecule_duplex_confirmed_passes() {
        let ref_ev = make_path_evidence(13, 0, 7, 6);
        // n_alt=1, duplex-confirmed. A duplex molecule contributes one read from
        // each strand, so min_simplex_fwd=1, min_simplex_rev=1.
        let alt_ev = PathEvidence {
            min_molecules: 1,
            mean_molecules: 1.0,
            min_duplex: 1,
            mean_duplex: 1.0,
            min_variant_specific_duplex: 1,
            mean_variant_specific_molecules: 1.0,
            min_simplex_fwd: 1,
            min_simplex_rev: 1,
            mean_error_prob: 0.001,
        };
        let cfg = CallerConfig::default();
        let call = call_variant("EGFR", &ref_ev, &alt_ev, b"G", b"A", &cfg);
        assert_eq!(call.filter, VariantFilter::Pass);
    }

    // Test 13: max_vaf filter rejects high-VAF calls.
    #[test]
    fn max_vaf_filters_germline() {
        let ref_ev = make_path_evidence(500, 0, 250, 250);
        // n_alt=500, VAF≈0.5 — likely germline.
        let alt_ev = make_path_evidence(500, 50, 250, 250);
        let cfg = CallerConfig {
            max_vaf: Some(0.35),
            ..CallerConfig::default()
        };
        let call = call_variant("EGFR", &ref_ev, &alt_ev, b"A", b"T", &cfg);
        assert_eq!(call.filter, VariantFilter::HighVaf);
    }

    // Test 14: large deletion (≥50 bp) is classified as LargeDeletion.
    #[test]
    fn classify_large_deletion() {
        let long_ref: Vec<u8> = b"ACGT".repeat(25).to_vec(); // 100 bp
        let short_alt: Vec<u8> = b"ACGT".repeat(10).to_vec(); // 40 bp — 60 bp deletion
        assert_eq!(
            classify_variant(&long_ref, &short_alt),
            VariantType::LargeDeletion
        );
    }

    // Test 15: small deletion (<50 bp) is still classified as Deletion.
    #[test]
    fn classify_small_deletion_not_sv() {
        let r: Vec<u8> = b"ACGT".repeat(15).to_vec(); // 60 bp
        let a: Vec<u8> = b"ACGT".repeat(10).to_vec(); // 40 bp — 20 bp deletion
        assert_eq!(classify_variant(&r, &a), VariantType::Deletion);
    }

    // Test 16: large insertion (≥50 bp) is classified as TandemDuplication.
    #[test]
    fn classify_tandem_duplication() {
        let short_ref: Vec<u8> = b"ACGT".repeat(10).to_vec(); // 40 bp
        let long_alt: Vec<u8> = b"ACGT".repeat(25).to_vec(); // 100 bp — 60 bp insertion
        assert_eq!(
            classify_variant(&short_ref, &long_alt),
            VariantType::TandemDuplication
        );
    }

    // Test 17: sequence that is its own reverse complement → not classified as Inversion.
    // Only genuine RC sequences (where alt == rc(ref)) get Inversion.
    #[test]
    fn classify_inversion() {
        // ref = AAATTT, rc = AAATTT — this is a palindrome, would be Inversion.
        // Use a non-palindrome: ref = ACGT, rc = ACGT — also palindrome.
        // Use ref = AACC, rc = GGTT — these differ.
        let ref_seq = b"AACCGGTT";
        // rc of AACCGGTT: complement T→A T→A G→C G→C C→G C→G A→T A→T = TTGGCCAA reversed
        // Actually let me compute manually:
        // AACCGGTT → complement: TTGGCCAA → reverse: AACCGGTT
        // That IS a palindrome. Let me use ACGTACGT.
        // ACGTACGT complement: TGCATGCA reverse: ACGTACGT — also palindrome.
        // ref = ACGTTGCA. rc = complement(ACGTTGCA) reversed = TGCAACGT reversed = TGCAACGT
        // complement: A→T C→G G→C T→A T→A G→C C→G A→T = TGCAACGT → reverse = TGCAACGT
        // Hmm let me use ref = AAACCC, rc:
        // complement: TTTGGG → reverse: GGGTTT ≠ AAACCC, so Inversion.
        let ref_seq = b"AAACCC";
        let alt_seq = b"GGGTTT"; // rc of AAACCC
        assert_eq!(classify_variant(ref_seq, alt_seq), VariantType::Inversion);
    }

    // Test 18: is_reverse_complement correctly identifies RC pairs.
    #[test]
    fn is_reverse_complement_basic() {
        // ACGT rc = ACGT (palindrome) — true.
        assert!(is_reverse_complement(b"ACGT", b"ACGT"));
        // AACC rc = GGTT — true.
        assert!(is_reverse_complement(b"AACC", b"GGTT"));
        // Different sequences, not RC.
        assert!(!is_reverse_complement(b"AACC", b"TTGG"));
        // Length 1 — returns false (too short).
        assert!(!is_reverse_complement(b"A", b"T"));
    }

    // Test 19: partial inversion in the central region with matching flanks.
    //
    // ref = AA [AAACCC] AA  (6 bp central, flanked by 2 bp each)
    // alt = AA [GGGTTT] AA  (central is rc of AAACCC)
    // Expected: Some(6) — 6 bp inverted segment.
    #[test]
    fn partial_inversion_len_central_segment() {
        // Flanks: AA ... AA.  Central 6 bp: AAACCC / GGGTTT.
        let ref_seq = b"AAAAACCCAA";
        let alt_seq = b"AAGGGTTTAA";
        assert_eq!(partial_inversion_len(ref_seq, alt_seq), Some(6));
    }

    // Test 20: partial inversion below SV_LENGTH_THRESHOLD is classified as MNV.
    //
    // A 6 bp inversion (< 50 bp) must not be promoted to VariantType::Inversion.
    #[test]
    fn partial_inversion_below_threshold_is_mnv() {
        // Same as test 19 — 6 bp inversion, below the 50 bp threshold.
        let ref_seq = b"AAAAACCCAA";
        let alt_seq = b"AAGGGTTTAA";
        assert_eq!(classify_variant(ref_seq, alt_seq), VariantType::Mnv);
    }

    // Test 21: full-path inversion still classified correctly after partial_inversion_len addition.
    //
    // Ensures the full-path RC check is not broken by the new code path.
    #[test]
    fn full_path_inversion_still_works() {
        // AAACCC rc = GGGTTT — full path RC, both checks should fire.
        let ref_seq = b"AAACCC";
        let alt_seq = b"GGGTTT";
        assert_eq!(classify_variant(ref_seq, alt_seq), VariantType::Inversion);
    }

    // Test 22: a single SNV inside an otherwise identical window → still SNV.
    //
    // partial_inversion_len must not confuse a single mismatch with an inversion.
    #[test]
    fn single_snv_not_misclassified_as_inversion() {
        let ref_seq = b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACGT"; // 103 bp
        let alt_seq = b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATGT"; // single A→T
        assert_eq!(classify_variant(ref_seq, alt_seq), VariantType::Snv);
    }

    // Test 23: SV type (LargeDeletion) uses sv_min_confidence, not min_confidence.
    //
    // A 100-bp deletion with 2 alt molecules gives confidence ~0.981.  Under
    // default SNV thresholds (min_confidence = 0.99) this would be LowConfidence.
    // With sv_min_confidence = 0.95 it should pass.
    #[test]
    fn sv_confidence_threshold_allows_low_confidence_del() {
        // 2 alt molecules out of 1002 total → confidence ≈ 0.981 < 0.99.
        let ref_ev = make_path_evidence(1000, 0, 500, 500);
        let alt_ev = make_path_evidence(2, 0, 1, 1);
        let long_ref: Vec<u8> = b"ACGT".repeat(25).to_vec(); // 100 bp
        let short_alt: Vec<u8> = b"ACGT".repeat(10).to_vec(); // 40 bp → 60 bp deletion
        let cfg = CallerConfig::default(); // sv_min_confidence = 0.95
        let call = call_variant("TP53", &ref_ev, &alt_ev, &long_ref, &short_alt, &cfg);
        // LargeDeletion → uses sv_min_confidence (0.95), passes at 0.981.
        assert_eq!(call.variant_type, VariantType::LargeDeletion);
        assert_eq!(call.filter, VariantFilter::Pass);
    }

    // Test 24: SV with confidence below sv_min_confidence is still LowConfidence.
    //
    // 1 alt molecule → confidence ≈ 0.83.  sv_min_confidence = 0.95 → LowConfidence.
    #[test]
    fn sv_confidence_threshold_rejects_very_low_confidence() {
        let ref_ev = make_path_evidence(1000, 0, 500, 500);
        let alt_ev = make_path_evidence(1, 0, 1, 0); // 1 mol, 0 duplex
        let long_ref: Vec<u8> = b"ACGT".repeat(25).to_vec();
        let short_alt: Vec<u8> = b"ACGT".repeat(10).to_vec();
        let cfg = CallerConfig {
            sv_min_alt_molecules: 1, // allow 1 mol for SVs
            ..CallerConfig::default()
        };
        let call = call_variant("TP53", &ref_ev, &alt_ev, &long_ref, &short_alt, &cfg);
        assert_eq!(call.variant_type, VariantType::LargeDeletion);
        // confidence ≈ 0.83 < sv_min_confidence 0.95 → LowConfidence.
        assert_eq!(call.filter, VariantFilter::LowConfidence);
    }

    // Test 25: SNV still uses min_confidence, not sv_min_confidence.
    //
    // 2 alt molecules → confidence ≈ 0.981 < 0.99, so LowConfidence for SNV
    // even though sv_min_confidence = 0.95 would allow it.
    #[test]
    fn snv_uses_standard_confidence_not_sv_threshold() {
        let ref_ev = make_path_evidence(1000, 0, 500, 500);
        let alt_ev = make_path_evidence(2, 0, 1, 1);
        let cfg = CallerConfig::default();
        let call = call_variant("KRAS", &ref_ev, &alt_ev, b"A", b"T", &cfg);
        assert_eq!(call.variant_type, VariantType::Snv);
        // confidence ≈ 0.981 < min_confidence 0.99 → LowConfidence.
        assert_eq!(call.filter, VariantFilter::LowConfidence);
    }

    // Test 26: INV with strand-biased alt reads passes when sv_strand_bias_threshold = 1.0.
    //
    // Inversion junction reads are structurally strand-biased: the de Bruijn
    // graph walks one orientation preferentially, so all alt simplex reads may
    // land on the forward strand. The default sv_strand_bias_threshold = 1.0
    // disables the filter for SV types, so the call must be PASS.
    #[test]
    fn inv_strand_biased_passes_with_default_sv_threshold() {
        // 60 bp inversion: alt is reverse complement of ref.
        let ref_seq: Vec<u8> = b"AAACCC".repeat(10).to_vec();
        let alt_seq: Vec<u8> = ref_seq
            .iter()
            .rev()
            .map(|&b| match b {
                b'A' => b'T',
                b'T' => b'A',
                b'C' => b'G',
                b'G' => b'C',
                other => other,
            })
            .collect();

        let ref_ev = make_path_evidence(990, 0, 495, 495);
        // Alt reads are all on the forward strand (structurally biased).
        let alt_ev = PathEvidence {
            min_molecules: 20,
            mean_molecules: 20.0,
            min_duplex: 5,
            mean_duplex: 5.0,
            min_variant_specific_duplex: 5,
            mean_variant_specific_molecules: 20.0,
            min_simplex_fwd: 20,
            min_simplex_rev: 0,
            mean_error_prob: 0.001,
        };

        // Default config: sv_strand_bias_threshold = 1.0 (disabled for SVs).
        let cfg = CallerConfig::default();
        let call = call_variant("BRCA1", &ref_ev, &alt_ev, &ref_seq, &alt_seq, &cfg);
        assert_eq!(call.variant_type, VariantType::Inversion);
        assert_eq!(
            call.filter,
            VariantFilter::Pass,
            "biased INV should pass when sv_strand_bias_threshold=1.0"
        );
    }

    // Test 27: same INV call is StrandBias-filtered when sv_strand_bias_threshold = 0.01.
    //
    // Setting sv_strand_bias_threshold to the SNV default re-enables the filter
    // for SV types, so an all-forward-strand alt must then be caught.
    #[test]
    fn inv_strand_biased_filtered_when_sv_threshold_set() {
        let ref_seq: Vec<u8> = b"AAACCC".repeat(10).to_vec();
        let alt_seq: Vec<u8> = ref_seq
            .iter()
            .rev()
            .map(|&b| match b {
                b'A' => b'T',
                b'T' => b'A',
                b'C' => b'G',
                b'G' => b'C',
                other => other,
            })
            .collect();

        let ref_ev = make_path_evidence(990, 0, 495, 495);
        let alt_ev = PathEvidence {
            min_molecules: 20,
            mean_molecules: 20.0,
            min_duplex: 5,
            mean_duplex: 5.0,
            min_variant_specific_duplex: 5,
            mean_variant_specific_molecules: 20.0,
            min_simplex_fwd: 20,
            min_simplex_rev: 0,
            mean_error_prob: 0.001,
        };

        // Lower sv_strand_bias_threshold to 0.01 to enable the filter for SVs.
        let cfg = CallerConfig {
            sv_strand_bias_threshold: 0.01,
            ..CallerConfig::default()
        };
        let call = call_variant("BRCA1", &ref_ev, &alt_ev, &ref_seq, &alt_seq, &cfg);
        assert_eq!(call.variant_type, VariantType::Inversion);
        assert_eq!(
            call.filter,
            VariantFilter::StrandBias,
            "biased INV should be filtered when sv_strand_bias_threshold=0.01"
        );
    }

    // Test 28: SNV with biased strands is still StrandBias-filtered.
    //
    // sv_strand_bias_threshold applies only to SV types. An SNV with all
    // alt reads on the forward strand must still be caught by strand_bias_threshold.
    #[test]
    fn snv_strand_biased_still_filtered() {
        let ref_ev = make_path_evidence(990, 0, 495, 495);
        let alt_ev = PathEvidence {
            min_molecules: 20,
            mean_molecules: 20.0,
            min_duplex: 5,
            mean_duplex: 5.0,
            min_variant_specific_duplex: 5,
            mean_variant_specific_molecules: 20.0,
            min_simplex_fwd: 20,
            min_simplex_rev: 0,
            mean_error_prob: 0.001,
        };

        // Default config: strand_bias_threshold = 0.01 for SNVs,
        // sv_strand_bias_threshold = 1.0 (disabled for SVs).
        let cfg = CallerConfig::default();
        let call = call_variant("TP53", &ref_ev, &alt_ev, b"A", b"T", &cfg);
        assert_eq!(call.variant_type, VariantType::Snv);
        assert_eq!(
            call.filter,
            VariantFilter::StrandBias,
            "biased SNV should still be filtered regardless of sv_strand_bias_threshold"
        );
    }

    // Test 29: SV path with high mean_variant_specific_molecules and low min_molecules
    //          uses the mean for NALT, not the minimum.
    //
    // A 100 bp deletion path might have min_molecules = 2 (one variant-specific
    // k-mer covered by only 2 molecules) but mean_variant_specific_molecules = 12.0
    // (the typical k-mer is covered by ~12 molecules at 1% VAF).  The call must
    // report n_molecules_alt = 12 (from the mean), not 2 (from the min).
    #[test]
    fn sv_call_uses_mean_variant_specific_not_min() {
        let long_ref: Vec<u8> = b"ACGT".repeat(25).to_vec(); // 100 bp
        let short_alt: Vec<u8> = b"ACGT".repeat(10).to_vec(); // 40 bp → 60 bp deletion
        assert_eq!(
            classify_variant(&long_ref, &short_alt),
            VariantType::LargeDeletion
        );

        let ref_ev = make_path_evidence(1000, 0, 500, 500);

        // min_molecules is low (bottleneck k-mer), but the mean is much higher.
        let alt_ev = PathEvidence {
            min_molecules: 2,
            mean_molecules: 12.0,
            min_duplex: 1,
            mean_duplex: 6.0,
            min_variant_specific_duplex: 1,
            mean_variant_specific_molecules: 12.0,
            min_simplex_fwd: 6,
            min_simplex_rev: 6,
            mean_error_prob: 0.001,
        };

        let cfg = CallerConfig::default();
        let call = call_variant("TP53", &ref_ev, &alt_ev, &long_ref, &short_alt, &cfg);
        assert_eq!(call.variant_type, VariantType::LargeDeletion);

        // n_molecules_alt should be 12 (mean), not 2 (min).
        assert_eq!(
            call.n_molecules_alt, 12,
            "SV call must use mean_variant_specific_molecules, got {}",
            call.n_molecules_alt
        );
        assert_eq!(call.filter, VariantFilter::Pass);
    }

    // Test 30: SNV call still uses min_molecules, not mean_variant_specific_molecules.
    //
    // For SNVs, mean_variant_specific_molecules is ignored. The call must use
    // min_molecules as before.
    #[test]
    fn snv_call_uses_min_molecules_unchanged() {
        let ref_ev = make_path_evidence(1000, 0, 500, 500);

        // min_molecules = 10, mean_variant_specific_molecules = 50.
        // For an SNV, k must be 10 (from min), not 50 (from mean).
        let alt_ev = PathEvidence {
            min_molecules: 10,
            mean_molecules: 50.0,
            min_duplex: 5,
            mean_duplex: 25.0,
            min_variant_specific_duplex: 5,
            mean_variant_specific_molecules: 50.0,
            min_simplex_fwd: 5,
            min_simplex_rev: 5,
            mean_error_prob: 0.001,
        };

        let cfg = CallerConfig::default();
        let call = call_variant("KRAS", &ref_ev, &alt_ev, b"A", b"T", &cfg);
        assert_eq!(call.variant_type, VariantType::Snv);

        // n_molecules_alt must equal min_molecules (10), not mean (50).
        assert_eq!(
            call.n_molecules_alt, 10,
            "SNV call must use min_molecules, got {}",
            call.n_molecules_alt
        );
        assert_eq!(call.filter, VariantFilter::Pass);
    }

    // Test 31 (BUG-004): estimate_vaf must not panic or wrap when k > m.
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
        assert!(lo >= 0.0 && lo <= 1.0, "CI lower out of range: {lo}");
        assert!(hi >= 0.0 && hi <= 1.0, "CI upper out of range: {hi}");
    }

    // Test 32 (FIX-004): call_variant with identical ref and alt must return a
    // filtered call, not a spurious Snv at 100% VAF.
    #[test]
    fn call_variant_identical_ref_alt_is_filtered() {
        let ref_ev = make_path_evidence(990, 0, 495, 495);
        let alt_ev = make_path_evidence(10, 5, 5, 5);
        let cfg = CallerConfig::default();
        let call = call_variant("KRAS", &ref_ev, &alt_ev, b"ACGT", b"ACGT", &cfg);
        assert_ne!(
            call.filter,
            VariantFilter::Pass,
            "identical ref/alt must not produce a PASS call"
        );
        assert_eq!(call.vaf, 0.0, "identical ref/alt should report VAF=0");
    }
}
