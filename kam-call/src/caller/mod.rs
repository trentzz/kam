//! Statistical variant calling using binomial/Beta-posterior models.
//!
//! Given molecule-level evidence from reference and alternate paths, this
//! module estimates VAF with 95% credible intervals, tests for strand bias,
//! and applies configurable filters to produce calibrated variant calls.
//!
//! The implementation is split across focused sub-modules:
//! - [`types`] — core types (`VariantCall`, `VariantFilter`, `VariantType`, `CallSource`)
//! - [`stats`] — statistical helpers (`estimate_vaf`, `strand_bias_test`, `compute_confidence`)
//! - [`classify`] — variant classification (`classify_variant` and sequence helpers)
//! - [`filter`] — filter assignment (`assign_filter`, `CallerConfig`)

mod classify;
mod filter;
mod stats;
mod types;

// Re-export everything to preserve the public API at `kam_call::caller::*`.
pub use classify::classify_variant;
pub use filter::{assign_filter, CallerConfig};
pub use stats::{compute_confidence, estimate_vaf, strand_bias_test};
pub use types::{CallSource, VariantCall, VariantFilter, VariantType, SV_LENGTH_THRESHOLD};

use kam_pathfind::score::PathEvidence;

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
            n_simplex_fwd_alt: 0,
            n_simplex_rev_alt: 0,
            n_duplex_ref: ref_evidence.min_duplex,
            n_simplex_ref: ref_evidence.min_simplex_fwd + ref_evidence.min_simplex_rev,
            mean_alt_error_prob: alt_evidence.mean_error_prob,
            min_variant_specific_duplex: 0,
            mean_variant_specific_molecules: 0.0,
            confidence: 0.0,
            strand_bias_p: 1.0,
            filter: VariantFilter::LowConfidence,
            ml_prob: None,
            call_source: CallSource::Called,
            rescue_min_alt_molecules: None,
            rescue_alt_duplex: None,
            rescue_approx_vaf: None,
            rescue_kmers_found: None,
            rescue_kmers_total: None,
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
        VariantType::LargeDeletion
            | VariantType::TandemDuplication
            | VariantType::Inversion
            | VariantType::Fusion
            | VariantType::InvDel
            | VariantType::NovelInsertion
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
        n_simplex_fwd_alt: alt_evidence.min_simplex_fwd,
        n_simplex_rev_alt: alt_evidence.min_simplex_rev,
        n_duplex_ref: ref_evidence.min_duplex,
        n_simplex_ref: ref_evidence.min_simplex_fwd + ref_evidence.min_simplex_rev,
        mean_alt_error_prob: alt_evidence.mean_error_prob,
        min_variant_specific_duplex: alt_evidence.min_variant_specific_duplex,
        mean_variant_specific_molecules: alt_evidence.mean_variant_specific_molecules,
        confidence,
        strand_bias_p,
        filter,
        ml_prob: None,
        call_source: CallSource::Called,
        rescue_min_alt_molecules: None,
        rescue_alt_duplex: None,
        rescue_approx_vaf: None,
        rescue_kmers_found: None,
        rescue_kmers_total: None,
    }
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
    //
    // No duplex alt molecules (min_variant_specific_duplex = 0), so the duplex
    // bypass does not engage and the biased call must be filtered.
    #[test]
    fn call_variant_strand_bias_flagged() {
        // Alt reads are all on the forward strand.
        let ref_ev = make_path_evidence(990, 0, 495, 495);
        let alt_ev = PathEvidence {
            min_molecules: 10,
            mean_molecules: 10.0,
            min_duplex: 0,
            mean_duplex: 0.0,
            min_variant_specific_duplex: 0, // no duplex → bypass does not apply
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

    // Test 48: integration test for InvDel — full pipeline from PathEvidence
    // through call_variant to VariantCall with correct type and PASS filter.
    #[test]
    fn integration_invdel_full_pipeline() {
        // Build an InvDel sequence pair (160 bp ref → 110 bp alt).
        let flank_l: Vec<u8> = b"AAAAAAAAAA".to_vec(); // 10 bp
        let del_region: Vec<u8> = b"TTTTTTTTTT".repeat(5).to_vec(); // 50 bp
        let inv_region: Vec<u8> = b"AAACCC".repeat(10).to_vec(); // 60 bp
        let flank_r: Vec<u8> = b"GGGGGGGGGG".repeat(4).to_vec(); // 40 bp

        let ref_seq: Vec<u8> = flank_l
            .iter()
            .chain(&del_region)
            .chain(&inv_region)
            .chain(&flank_r)
            .copied()
            .collect();

        let rc_inv: Vec<u8> = inv_region
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
        let alt_seq: Vec<u8> = flank_l
            .iter()
            .chain(&rc_inv)
            .chain(&flank_r)
            .copied()
            .collect();

        // Verify the classification.
        assert_eq!(classify_variant(&ref_seq, &alt_seq), VariantType::InvDel);

        // Construct PathEvidence with SV-typical characteristics:
        // low min_molecules (long path bottleneck) but good mean_variant_specific.
        let ref_ev = PathEvidence {
            min_molecules: 1000,
            mean_molecules: 1000.0,
            min_duplex: 0,
            mean_duplex: 0.0,
            min_variant_specific_duplex: 0,
            mean_variant_specific_molecules: 1000.0,
            min_simplex_fwd: 500,
            min_simplex_rev: 500,
            mean_error_prob: 0.001,
        };
        let alt_ev = PathEvidence {
            min_molecules: 2, // bottlenecked at a rare k-mer
            mean_molecules: 15.0,
            min_duplex: 1,
            mean_duplex: 7.0,
            min_variant_specific_duplex: 1,
            mean_variant_specific_molecules: 15.0, // 1.5% VAF estimate
            min_simplex_fwd: 8,
            min_simplex_rev: 7,
            mean_error_prob: 0.001,
        };

        let cfg = CallerConfig::default();
        let call = call_variant("BRCA2", &ref_ev, &alt_ev, &ref_seq, &alt_seq, &cfg);

        assert_eq!(
            call.variant_type,
            VariantType::InvDel,
            "type must be InvDel"
        );
        assert_eq!(call.filter, VariantFilter::Pass, "InvDel must pass filters");
        // n_molecules_alt for SV types uses mean_variant_specific_molecules.
        assert_eq!(
            call.n_molecules_alt, 15,
            "n_molecules_alt must use mean_variant_specific_molecules"
        );
        assert!(
            call.vaf > 0.0 && call.vaf < 0.05,
            "VAF must be small (around 1.5%), got {}",
            call.vaf
        );
    }

    // Test 49: integration test for NovelInsertion — full pipeline.
    #[test]
    fn integration_novel_insertion_full_pipeline() {
        let ref_seq: Vec<u8> = b"ACGT".repeat(25).to_vec(); // 100 bp
                                                            // Insert 60 bp of all-C at position 50 — not present in ACGT repeats.
        let novel_insert: Vec<u8> = vec![b'C'; 60];
        let alt_seq: Vec<u8> = ref_seq[..50]
            .iter()
            .chain(&novel_insert)
            .chain(&ref_seq[50..])
            .copied()
            .collect();

        assert_eq!(
            classify_variant(&ref_seq, &alt_seq),
            VariantType::NovelInsertion
        );

        let ref_ev = PathEvidence {
            min_molecules: 800,
            mean_molecules: 800.0,
            min_duplex: 0,
            mean_duplex: 0.0,
            min_variant_specific_duplex: 0,
            mean_variant_specific_molecules: 800.0,
            min_simplex_fwd: 400,
            min_simplex_rev: 400,
            mean_error_prob: 0.001,
        };
        let alt_ev = PathEvidence {
            min_molecules: 1,
            mean_molecules: 10.0,
            min_duplex: 0,
            mean_duplex: 5.0,
            min_variant_specific_duplex: 0,
            mean_variant_specific_molecules: 10.0,
            min_simplex_fwd: 5,
            min_simplex_rev: 5,
            mean_error_prob: 0.001,
        };

        let cfg = CallerConfig::default();
        let call = call_variant("ALK", &ref_ev, &alt_ev, &ref_seq, &alt_seq, &cfg);

        assert_eq!(
            call.variant_type,
            VariantType::NovelInsertion,
            "type must be NovelInsertion"
        );
        assert_eq!(
            call.filter,
            VariantFilter::Pass,
            "NovelInsertion must pass filters"
        );
    }
}
