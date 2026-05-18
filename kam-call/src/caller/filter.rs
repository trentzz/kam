//! Filter assignment and caller configuration.

use super::types::{VariantFilter, VariantType};

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

/// Assign a filter label given statistical evidence and thresholds.
///
/// SV types (`LargeDeletion`, `TandemDuplication`, `Inversion`) use
/// `sv_min_confidence` and `sv_min_alt_molecules` instead of the SNV/indel
/// defaults. A large structural event requires fewer supporting molecules to
/// reach the same confidence as a single-base change, because background
/// sequencing errors cannot produce a 50–200 bp structural signature.
pub fn assign_filter(
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
        VariantType::LargeDeletion
            | VariantType::TandemDuplication
            | VariantType::Inversion
            | VariantType::Fusion
            | VariantType::InvDel
            | VariantType::NovelInsertion
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
    // Duplex molecules are sequenced on both strands by definition, so a call
    // supported by ≥2 duplex alt molecules cannot have genuine strand bias.
    // Skip the filter in that case to avoid penalising real variants whose
    // apparent bias is driven by reference-strand asymmetry rather than artefact.
    let duplex_bypass = n_duplex_alt >= 2;
    if !duplex_bypass
        && eff_strand_bias_threshold < 1.0
        && strand_bias_p < eff_strand_bias_threshold
    {
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
    use super::super::stats::compute_confidence;
    use super::*;

    // Test: CallerConfig default values match documented expectations.
    // The defaults are public API and must not drift silently.
    #[test]
    fn caller_config_default_values() {
        let cfg = CallerConfig::default();
        assert_eq!(cfg.min_confidence, 0.99, "min_confidence");
        assert_eq!(cfg.strand_bias_threshold, 0.01, "strand_bias_threshold");
        assert_eq!(cfg.min_alt_molecules, 2, "min_alt_molecules");
        assert_eq!(cfg.min_alt_duplex, 0, "min_alt_duplex");
        assert_eq!(cfg.min_alt_duplex_for_single, Some(1), "duplex_for_single");
        assert_eq!(cfg.max_vaf, None, "max_vaf");
        assert_eq!(cfg.sv_min_confidence, 0.95, "sv_min_confidence");
        assert_eq!(cfg.sv_min_alt_molecules, 1, "sv_min_alt_molecules");
        assert_eq!(
            cfg.sv_strand_bias_threshold, 1.0,
            "sv_strand_bias_threshold"
        );
    }

    // Test 35: assign_filter applies SV thresholds to Fusion, InvDel, and NovelInsertion.
    #[test]
    fn new_sv_types_use_sv_thresholds_in_assign_filter() {
        // 2 alt molecules → confidence ≈ 0.981, which passes sv_min_confidence = 0.95
        // but fails min_confidence = 0.99.
        let confidence = compute_confidence(2, 1002, 1e-4);
        assert!(
            confidence < 0.99 && confidence > 0.95,
            "confidence {confidence} must be between 0.95 and 0.99 for this test to be meaningful"
        );
        let cfg = CallerConfig::default();

        for vt in [
            VariantType::Fusion,
            VariantType::InvDel,
            VariantType::NovelInsertion,
        ] {
            let filter = assign_filter(
                confidence, 0.5,   // no strand bias
                2,     // n_alt = 2, above sv_min_alt_molecules = 1
                0,     // n_duplex_alt = 0
                0.002, // vaf
                vt, &cfg,
            );
            assert_eq!(
                filter,
                VariantFilter::Pass,
                "{vt} should PASS at confidence {confidence:.4} with SV thresholds"
            );
        }

        // Confirm the same confidence is LowConfidence for a non-SV type.
        let snv_filter = assign_filter(confidence, 0.5, 2, 0, 0.002, VariantType::Snv, &cfg);
        assert_eq!(
            snv_filter,
            VariantFilter::LowConfidence,
            "SNV must use min_confidence (0.99), not sv_min_confidence"
        );
    }

    // Test 45: a call with confidence 0.96 passes for all three new SV types
    // (sv_min_confidence = 0.95) but would fail for an SNV (min_confidence = 0.99).
    #[test]
    fn sv_confidence_96_passes_new_types_fails_snv() {
        let confidence = 0.96_f64;
        let cfg = CallerConfig::default();

        for vt in [
            VariantType::Fusion,
            VariantType::InvDel,
            VariantType::NovelInsertion,
        ] {
            let filter = assign_filter(
                confidence, 0.5,   // strand bias p — balanced, no flag
                2,     // n_alt >= sv_min_alt_molecules (1)
                0,     // n_duplex_alt
                0.002, // vaf
                vt, &cfg,
            );
            assert_eq!(
                filter,
                VariantFilter::Pass,
                "{vt} at confidence 0.96 must PASS with sv_min_confidence=0.95"
            );
        }

        // Same confidence must be LowConfidence for SNV.
        let snv_filter = assign_filter(confidence, 0.5, 2, 0, 0.002, VariantType::Snv, &cfg);
        assert_eq!(
            snv_filter,
            VariantFilter::LowConfidence,
            "SNV at confidence 0.96 must fail min_confidence=0.99"
        );
    }

    // Test 46: all three new SV types use sv_min_alt_molecules (1), not
    // min_alt_molecules (2). A single alt molecule without duplex must pass.
    #[test]
    fn new_sv_types_pass_with_one_alt_molecule() {
        let confidence = 0.97_f64;
        let cfg = CallerConfig {
            // Disable single-molecule duplex bypass to test the sv threshold directly.
            min_alt_duplex_for_single: None,
            ..CallerConfig::default()
        };

        for vt in [
            VariantType::Fusion,
            VariantType::InvDel,
            VariantType::NovelInsertion,
        ] {
            let filter = assign_filter(
                confidence, 0.5, // no strand bias
                1,   // n_alt = 1 — below SNV threshold (2) but at SV threshold (1)
                0,   // no duplex
                0.001, vt, &cfg,
            );
            assert_eq!(
                filter,
                VariantFilter::Pass,
                "{vt} with n_alt=1 must PASS using sv_min_alt_molecules=1"
            );
        }

        // SNV with n_alt=1 and no duplex must be LowConfidence.
        let snv_filter = assign_filter(confidence, 0.5, 1, 0, 0.001, VariantType::Snv, &cfg);
        assert_eq!(
            snv_filter,
            VariantFilter::LowConfidence,
            "SNV with n_alt=1 and no duplex must fail min_alt_molecules=2"
        );
    }

    // Test 47: all three new SV types respect sv_strand_bias_threshold = 1.0
    // (disabled by default). A strand-biased alt passes for SV but fails for SNV.
    #[test]
    fn new_sv_types_ignore_strand_bias_by_default() {
        // p_strand = 0.001 — very biased; below SNV threshold (0.01).
        let p_strand = 0.001_f64;
        let confidence = 0.97_f64;
        let cfg = CallerConfig::default(); // sv_strand_bias_threshold = 1.0 (disabled)

        for vt in [
            VariantType::Fusion,
            VariantType::InvDel,
            VariantType::NovelInsertion,
        ] {
            let filter = assign_filter(confidence, p_strand, 2, 0, 0.002, vt, &cfg);
            assert_eq!(
                filter,
                VariantFilter::Pass,
                "{vt} must ignore strand bias when sv_strand_bias_threshold=1.0"
            );
        }

        // SNV with same biased strand must be filtered.
        let snv_filter = assign_filter(confidence, p_strand, 2, 0, 0.002, VariantType::Snv, &cfg);
        assert_eq!(
            snv_filter,
            VariantFilter::StrandBias,
            "SNV must still be filtered by strand bias"
        );
    }

    // Test 59: strand-biased deletion passes when n_duplex_alt >= 2 (duplex bypass).
    #[test]
    fn strand_bias_bypassed_when_duplex_alt_sufficient() {
        let cfg = CallerConfig {
            min_confidence: 0.80,
            ..CallerConfig::default()
        };
        // strand_bias_p = 0.005, below default threshold of 0.01 → normally StrandBias.
        // n_duplex_alt = 2 → bypass must engage and return Pass.
        let filter = assign_filter(0.95, 0.005, 12, 2, 0.03, VariantType::Deletion, &cfg);
        assert_eq!(
            filter,
            VariantFilter::Pass,
            "duplex-confirmed call must not be filtered for strand bias"
        );
    }

    // Test 60: strand-biased deletion is still filtered when n_duplex_alt < 2.
    #[test]
    fn strand_bias_still_filtered_when_duplex_alt_insufficient() {
        let cfg = CallerConfig {
            min_confidence: 0.80,
            ..CallerConfig::default()
        };
        // strand_bias_p = 0.005, below default threshold of 0.01 → StrandBias.
        // n_duplex_alt = 1 → bypass does not engage.
        let filter = assign_filter(0.95, 0.005, 12, 1, 0.03, VariantType::Deletion, &cfg);
        assert_eq!(
            filter,
            VariantFilter::StrandBias,
            "call with only 1 duplex alt must still be filtered for strand bias"
        );
    }
}
