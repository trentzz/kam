//! Optional ML scoring via ONNX model.
//!
//! When an ONNX model and companion `model_meta.json` are provided,
//! each [`VariantCall`] is scored by a gradient-boosted classifier
//! trained on simulated training data. The resulting probability
//! is stored in [`VariantCall::ml_prob`].
//!
//! # Example
//!
//! ```no_run
//! use std::path::Path;
//! use kam_call::ml::MlScorer;
//!
//! let scorer = MlScorer::load(
//!     Path::new("model.onnx"),
//!     Path::new("model_meta.json"),
//! ).expect("load model");
//! ```

#[cfg(feature = "ml")]
use ndarray::Array2;
#[cfg(feature = "ml")]
use ort::{
    session::{builder::GraphOptimizationLevel, Session},
    value::TensorRef,
};

use std::collections::HashMap;
use std::path::Path;

use crate::caller::VariantCall;
#[cfg(feature = "ml")]
use crate::caller::VariantType;

/// Metadata loaded alongside the ONNX model.
#[derive(Debug, serde::Deserialize)]
pub struct ModelMeta {
    /// Feature names in the exact order the model expects.
    pub feature_names: Vec<String>,
    /// Probability threshold above which a call is labelled ML_PASS.
    pub ml_pass_threshold: f64,
    /// Mapping from variant type string to integer encoding.
    pub variant_class_map: HashMap<String, i32>,
}

/// Holds a loaded ONNX session and associated metadata.
pub struct MlScorer {
    #[cfg(feature = "ml")]
    session: Session,
    pub meta: ModelMeta,
}

impl MlScorer {
    /// Load an ONNX model and its companion metadata file.
    ///
    /// # Errors
    ///
    /// Returns an error if either file cannot be read or parsed.
    pub fn load(model_path: &Path, meta_path: &Path) -> Result<Self, Box<dyn std::error::Error>> {
        let meta_bytes = std::fs::read(meta_path)?;

        #[cfg(feature = "ml")]
        {
            let meta: ModelMeta = serde_json::from_slice(&meta_bytes)?;
            let session = Session::builder()?
                .with_optimization_level(GraphOptimizationLevel::Level3)?
                .with_intra_threads(1)?
                .commit_from_file(model_path)?;
            Ok(Self { session, meta })
        }

        #[cfg(not(feature = "ml"))]
        {
            let _ = (model_path, meta_bytes);
            Err("kam-call was built without the 'ml' feature — recompile with --features ml".into())
        }
    }

    /// Load an ONNX model and metadata from in-memory bytes (e.g. from `include_bytes!`).
    ///
    /// # Errors
    ///
    /// Returns an error if the bytes cannot be parsed.
    pub fn from_bytes(
        model_bytes: &[u8],
        meta_bytes: &[u8],
    ) -> Result<Self, Box<dyn std::error::Error>> {
        #[cfg(feature = "ml")]
        {
            let meta: ModelMeta = serde_json::from_slice(meta_bytes)?;
            let session = Session::builder()?
                .with_optimization_level(GraphOptimizationLevel::Level3)?
                .with_intra_threads(1)?
                .commit_from_memory(model_bytes)?;
            Ok(Self { session, meta })
        }

        #[cfg(not(feature = "ml"))]
        {
            let _ = (model_bytes, meta_bytes);
            Err("kam-call was built without the 'ml' feature — recompile with --features ml".into())
        }
    }

    /// Score a single [`VariantCall`] and return the ML probability (class 1).
    ///
    /// Returns `None` if feature extraction or inference fails.
    #[cfg(feature = "ml")]
    pub fn score(&mut self, call: &VariantCall) -> Option<f32> {
        let features = self.extract_features(call);
        let n = self.meta.feature_names.len();
        if features.len() != n {
            return None;
        }

        let array = Array2::from_shape_vec((1, n), features).ok()?;
        let tensor = TensorRef::from_array_view(array.view()).ok()?;
        let outputs = self.session.run(ort::inputs![tensor]).ok()?;

        // LightGBM ONNX output: "label" (i64) and "probabilities" (shape [N, 2]).
        // Index 1 of the second axis is P(positive).
        let probs = outputs
            .get("probabilities")
            .or_else(|| outputs.get("output_probability"))?;
        let (_shape, data) = probs.try_extract_tensor::<f32>().ok()?;
        // data has layout [n_samples * n_classes]; index 1 is P(class 1) for sample 0.
        data.get(1).copied()
    }

    /// Score a single [`VariantCall`]. Returns `None` when the ml feature is not compiled in.
    #[cfg(not(feature = "ml"))]
    pub fn score(&mut self, _call: &VariantCall) -> Option<f32> {
        None
    }

    /// Build the feature vector for a call in the order `meta.feature_names` specifies.
    #[cfg(feature = "ml")]
    fn extract_features(&self, call: &VariantCall) -> Vec<f32> {
        // ── Original 33 features ──────────────────────────────────────────────
        let vaf = call.vaf as f32;
        let nref = call.n_molecules_ref as f32;
        let nalt = call.n_molecules_alt as f32;
        let ndupalt = call.n_duplex_alt as f32;
        let nsimalt = call.n_simplex_alt as f32;
        let sbp = call.strand_bias_p as f32;
        let conf = call.confidence as f32;
        let ref_len = call.ref_sequence.len() as f32;
        let alt_len = call.alt_sequence.len() as f32;

        let duplex_frac = ndupalt / (nalt + 1e-9);
        let has_duplex = if ndupalt > 0.0 { 1.0_f32 } else { 0.0_f32 };
        let ci_width = (call.vaf_ci_high - call.vaf_ci_low) as f32;
        let alt_depth = nref + nalt;

        let log_nalt = (nalt + 1.0).ln();
        let log_nref = (nref + 1.0).ln();
        let log_alt_depth = (alt_depth + 1.0).ln();
        let log_vaf = (vaf + 1e-6).ln();

        let vaf_times_conf = vaf * conf;
        let vaf_times_nalt = vaf * nalt;
        let nalt_over_conf = nalt / (conf + 1e-9);
        let ci_width_rel = ci_width / (vaf + 1e-9);
        let snr = nalt / (nref + 1.0);

        let conf_sq = conf * conf;
        let nalt_sq = nalt * nalt;
        let vaf_sq = vaf * vaf;

        let ref_alt_len_ratio = ref_len / (alt_len + 1e-9_f32);
        let indel_size = (ref_len - alt_len).abs();

        let duplex_enrichment = ndupalt / (vaf * alt_depth + 1e-9);
        let simplex_only_frac = nsimalt / (nalt + 1e-9);

        let conf_above_99 = if conf >= 0.99 { 1.0_f32 } else { 0.0_f32 };
        let conf_above_999 = if conf >= 0.999 { 1.0_f32 } else { 0.0_f32 };
        let sbp_above_05 = if sbp >= 0.05 { 1.0_f32 } else { 0.0_f32 };

        let variant_class_enc = *self
            .meta
            .variant_class_map
            .get(variant_type_str(call.variant_type))
            .unwrap_or(&0) as f32;

        // ── Category B: existing VariantCall fields not previously exposed ────
        let n_simplex_fwd_alt = call.n_simplex_fwd_alt as f32;
        let n_simplex_rev_alt = call.n_simplex_rev_alt as f32;
        let n_duplex_ref = call.n_duplex_ref as f32;
        let n_simplex_ref = call.n_simplex_ref as f32;
        let mean_alt_error_prob = call.mean_alt_error_prob;
        let min_variant_specific_duplex = call.min_variant_specific_duplex as f32;
        let mean_variant_specific_molecules = call.mean_variant_specific_molecules;

        // ── Category C: derived strand/duplex features ────────────────────────
        let strand_asymmetry_alt = (n_simplex_fwd_alt - n_simplex_rev_alt)
            / (n_simplex_fwd_alt + n_simplex_rev_alt + 1e-9);
        let duplex_vaf = ndupalt / (ndupalt + n_duplex_ref + 1e-9);
        let simplex_vaf = nsimalt / (nsimalt + n_simplex_ref + 1e-9);
        let duplex_simplex_vaf_delta = duplex_vaf - simplex_vaf;

        // ── Category A: sequence-context features ─────────────────────────────
        let subst_type = snv_subst_type(&call.ref_sequence, &call.alt_sequence, call.variant_type);
        let trinuc_context =
            snv_trinuc_context(&call.ref_sequence, &call.alt_sequence, call.variant_type);
        let is_cpg = snv_is_cpg(&call.ref_sequence, &call.alt_sequence, call.variant_type);
        let gc_content_ref = ref_gc_content(&call.ref_sequence);
        let homopolymer_run =
            ref_homopolymer_run(&call.ref_sequence, &call.alt_sequence, call.variant_type);

        let lookup: HashMap<&str, f32> = [
            // original 33
            ("vaf", vaf),
            ("nref", nref),
            ("nalt", nalt),
            ("ndupalt", ndupalt),
            ("nsimalt", nsimalt),
            ("sbp", sbp),
            ("conf", conf),
            ("ref_len", ref_len),
            ("alt_len", alt_len),
            ("duplex_frac", duplex_frac),
            ("has_duplex", has_duplex),
            ("ci_width", ci_width),
            ("alt_depth", alt_depth),
            ("log_nalt", log_nalt),
            ("log_nref", log_nref),
            ("log_alt_depth", log_alt_depth),
            ("log_vaf", log_vaf),
            ("vaf_times_conf", vaf_times_conf),
            ("vaf_times_nalt", vaf_times_nalt),
            ("nalt_over_conf", nalt_over_conf),
            ("ci_width_rel", ci_width_rel),
            ("snr", snr),
            ("conf_sq", conf_sq),
            ("nalt_sq", nalt_sq),
            ("vaf_sq", vaf_sq),
            ("ref_alt_len_ratio", ref_alt_len_ratio),
            ("indel_size", indel_size),
            ("duplex_enrichment", duplex_enrichment),
            ("simplex_only_frac", simplex_only_frac),
            ("conf_above_99", conf_above_99),
            ("conf_above_999", conf_above_999),
            ("sbp_above_05", sbp_above_05),
            ("variant_class_enc", variant_class_enc),
            // category B
            ("n_simplex_fwd_alt", n_simplex_fwd_alt),
            ("n_simplex_rev_alt", n_simplex_rev_alt),
            ("n_duplex_ref", n_duplex_ref),
            ("n_simplex_ref", n_simplex_ref),
            ("mean_alt_error_prob", mean_alt_error_prob),
            ("min_variant_specific_duplex", min_variant_specific_duplex),
            (
                "mean_variant_specific_molecules",
                mean_variant_specific_molecules,
            ),
            // category C
            ("strand_asymmetry_alt", strand_asymmetry_alt),
            ("duplex_vaf", duplex_vaf),
            ("simplex_vaf", simplex_vaf),
            ("duplex_simplex_vaf_delta", duplex_simplex_vaf_delta),
            // category A
            ("subst_type", subst_type),
            ("trinuc_context", trinuc_context),
            ("is_cpg", is_cpg),
            ("gc_content_ref", gc_content_ref),
            ("homopolymer_run", homopolymer_run),
            ("dust_score", compute_dust_score(&call.ref_sequence, 64)),
            ("repeat_fraction", compute_repeat_fraction(&call.ref_sequence)),
        ]
        .into_iter()
        .collect();

        self.meta
            .feature_names
            .iter()
            .map(|name| *lookup.get(name.as_str()).unwrap_or(&0.0))
            .collect()
    }
}

#[cfg(feature = "ml")]
fn variant_type_str(vt: VariantType) -> &'static str {
    match vt {
        VariantType::Snv => "SNV",
        VariantType::Insertion => "Insertion",
        VariantType::Deletion => "Deletion",
        VariantType::Mnv => "MNV",
        VariantType::Complex => "Complex",
        VariantType::LargeDeletion => "LargeDeletion",
        VariantType::TandemDuplication => "TandemDuplication",
        VariantType::Inversion => "Inversion",
        VariantType::Fusion => "Fusion",
        VariantType::InvDel => "InvDel",
        VariantType::NovelInsertion => "NovelInsertion",
    }
}

/// Encode the SNV substitution class (0–11) or return 12 for non-SNV/unrecognised.
///
/// Uses 12 canonical classes:
/// C>A=0, C>G=1, C>T=2, T>A=3, T>C=4, T>G=5, G>T=6, G>C=7, G>A=8, A>T=9, A>G=10, A>C=11.
#[cfg(feature = "ml")]
fn snv_subst_type(ref_seq: &[u8], alt_seq: &[u8], vt: VariantType) -> f32 {
    if vt != VariantType::Snv {
        return 12.0;
    }
    let diff = ref_seq.iter().zip(alt_seq.iter()).position(|(r, a)| r != a);
    let i = match diff {
        Some(i) => i,
        None => return 12.0,
    };
    match (ref_seq[i], alt_seq[i]) {
        (b'C', b'A') => 0.0,
        (b'C', b'G') => 1.0,
        (b'C', b'T') => 2.0,
        (b'T', b'A') => 3.0,
        (b'T', b'C') => 4.0,
        (b'T', b'G') => 5.0,
        (b'G', b'T') => 6.0,
        (b'G', b'C') => 7.0,
        (b'G', b'A') => 8.0,
        (b'A', b'T') => 9.0,
        (b'A', b'G') => 10.0,
        (b'A', b'C') => 11.0,
        _ => 12.0,
    }
}

/// Encode the trinucleotide context centred on the SNV site (0–63), or 64 for non-SNV/edge.
///
/// Uses base-4 encoding: A=0, C=1, G=2, T=3.
/// Value = ref[i-1]*16 + ref[i]*4 + ref[i+1].
#[cfg(feature = "ml")]
fn snv_trinuc_context(ref_seq: &[u8], alt_seq: &[u8], vt: VariantType) -> f32 {
    if vt != VariantType::Snv {
        return 64.0;
    }
    let diff = ref_seq.iter().zip(alt_seq.iter()).position(|(r, a)| r != a);
    let i = match diff {
        Some(i) if i >= 1 && i + 1 < ref_seq.len() => i,
        _ => return 64.0,
    };
    let encode = |b: u8| -> Option<u8> {
        match b {
            b'A' | b'a' => Some(0),
            b'C' | b'c' => Some(1),
            b'G' | b'g' => Some(2),
            b'T' | b't' => Some(3),
            _ => None,
        }
    };
    match (
        encode(ref_seq[i - 1]),
        encode(ref_seq[i]),
        encode(ref_seq[i + 1]),
    ) {
        (Some(l), Some(c), Some(r)) => (l as f32) * 16.0 + (c as f32) * 4.0 + (r as f32),
        _ => 64.0,
    }
}

/// Return 1.0 if the SNV is a C>T or G>A change at a CpG dinucleotide, else 0.0.
///
/// C>T is CpG-context when the base immediately 3' of the C is G.
/// G>A is CpG-context when the base immediately 5' of the G is C.
#[cfg(feature = "ml")]
fn snv_is_cpg(ref_seq: &[u8], alt_seq: &[u8], vt: VariantType) -> f32 {
    if vt != VariantType::Snv {
        return 0.0;
    }
    let diff = ref_seq.iter().zip(alt_seq.iter()).position(|(r, a)| r != a);
    let i = match diff {
        Some(i) => i,
        None => return 0.0,
    };
    let is_ct = ref_seq[i] == b'C' && alt_seq[i] == b'T';
    let is_ga = ref_seq[i] == b'G' && alt_seq[i] == b'A';
    if is_ct && i + 1 < ref_seq.len() && ref_seq[i + 1] == b'G' {
        return 1.0;
    }
    if is_ga && i >= 1 && ref_seq[i - 1] == b'C' {
        return 1.0;
    }
    0.0
}

/// Return the GC fraction of the reference sequence (0.0 for empty input).
#[cfg(feature = "ml")]
fn ref_gc_content(ref_seq: &[u8]) -> f32 {
    if ref_seq.is_empty() {
        return 0.0;
    }
    let gc = ref_seq
        .iter()
        .filter(|&&b| b == b'G' || b == b'g' || b == b'C' || b == b'c')
        .count();
    gc as f32 / ref_seq.len() as f32
}

/// Return the length of the longest homopolymer run adjacent to the variant site.
///
/// Finds the first differing position and scans left and right in `ref_seq`
/// counting runs of the same base. Returns the longer of the two runs.
/// Returns 0.0 for non-SNV or when the variant position is at an edge.
#[cfg(feature = "ml")]
fn ref_homopolymer_run(ref_seq: &[u8], alt_seq: &[u8], vt: VariantType) -> f32 {
    if vt != VariantType::Snv {
        return 0.0;
    }
    let diff = ref_seq.iter().zip(alt_seq.iter()).position(|(r, a)| r != a);
    let i = match diff {
        Some(i) => i,
        None => return 0.0,
    };
    // Scan left from i-1.
    let left_run = if i > 0 {
        let base = ref_seq[i - 1];
        ref_seq[..i]
            .iter()
            .rev()
            .take_while(|&&b| b == base)
            .count()
    } else {
        0
    };
    // Scan right from i+1.
    let right_run = if i + 1 < ref_seq.len() {
        let base = ref_seq[i + 1];
        ref_seq[i + 1..].iter().take_while(|&&b| b == base).count()
    } else {
        0
    };
    left_run.max(right_run) as f32
}

/// Compute a DUST-like sequence complexity score for the reference sequence.
/// Higher values indicate more repetitive, low-complexity sequence.
/// A score above 10 suggests borderline low complexity; above 30 is high repeat.
///
/// # Examples
/// ```ignore
/// // Poly-A sequence scores high (low complexity)
/// let score = compute_dust_score(b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", 64);
/// assert!(score > 10.0);
/// // Random-looking sequence scores low
/// let score = compute_dust_score(b"ACGTGCTAGCTAGCATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG", 64);
/// assert!(score < 10.0);
/// ```
#[cfg(feature = "ml")]
fn compute_dust_score(seq: &[u8], window: usize) -> f32 {
    if seq.len() < 3 {
        return 0.0;
    }
    let mut max_score: f32 = 0.0;
    let window = window.min(seq.len());
    let n_windows = seq.len() - window + 1;
    for start in 0..n_windows {
        let end = start + window;
        let w = &seq[start..end];
        // Count all trinucleotides in this window.
        let mut counts = [0u32; 64];
        for i in 0..w.len().saturating_sub(2) {
            let a = base_index(w[i]);
            let b = base_index(w[i + 1]);
            let c = base_index(w[i + 2]);
            let idx = a * 16 + b * 4 + c;
            counts[idx] += 1;
        }
        let denom = (w.len().saturating_sub(2)).max(1) as f32;
        let score: f32 = counts
            .iter()
            .map(|&c| (c * c.saturating_sub(1) / 2) as f32)
            .sum::<f32>()
            / denom;
        if score > max_score {
            max_score = score;
        }
    }
    max_score
}

/// Map a nucleotide byte to 0..=3 for trinucleotide indexing.
/// Non-ACGT bases map to 0 (A).
#[cfg(feature = "ml")]
fn base_index(b: u8) -> usize {
    match b.to_ascii_uppercase() {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        b'T' => 3,
        _ => 0,
    }
}

/// Compute the fraction of bases in the reference sequence that fall within
/// homopolymer runs (≥3 identical bases) or dinucleotide repeat runs (≥6 bases).
///
/// # Examples
/// ```ignore
/// // Poly-A has high repeat fraction
/// let f = compute_repeat_fraction(b"AAAAAACGT");
/// assert!(f > 0.5);
/// // Short sequence with no repeats
/// let f = compute_repeat_fraction(b"ACGT");
/// assert_eq!(f, 0.0);
/// ```
#[cfg(feature = "ml")]
fn compute_repeat_fraction(seq: &[u8]) -> f32 {
    if seq.is_empty() {
        return 0.0;
    }
    let n = seq.len();
    let mut in_repeat = vec![false; n];

    // Homopolymer runs >= 3.
    let mut i = 0;
    while i < n {
        let mut j = i + 1;
        while j < n && seq[j] == seq[i] {
            j += 1;
        }
        if j - i >= 3 {
            in_repeat[i..j].iter_mut().for_each(|b| *b = true);
        }
        i = j;
    }

    // Dinucleotide repeats >= 3 units (>= 6 bases).
    if n >= 6 {
        for start in 0..n.saturating_sub(5) {
            let di = &seq[start..start + 2];
            let mut end = start + 2;
            while end + 1 < n && seq[end] == di[0] && seq[end + 1] == di[1] {
                end += 2;
            }
            if end - start >= 6 {
                in_repeat[start..end].iter_mut().for_each(|b| *b = true);
            }
        }
    }

    let repeat_count = in_repeat.iter().filter(|&&b| b).count();
    repeat_count as f32 / n as f32
}

#[cfg(all(test, feature = "ml"))]
mod tests {
    use super::*;

    #[test]
    fn dust_score_poly_a_is_high() {
        let seq = b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
        assert!(
            compute_dust_score(seq, 64) > 10.0,
            "poly-A should have high DUST score"
        );
    }

    #[test]
    fn dust_score_complex_sequence_is_low() {
        let seq = b"ACGTGCTAGCTAGCATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG";
        assert!(
            compute_dust_score(seq, 64) < 10.0,
            "complex sequence should have low DUST score"
        );
    }

    #[test]
    fn dust_score_empty_sequence_returns_zero() {
        assert_eq!(compute_dust_score(b"", 64), 0.0);
    }

    #[test]
    fn repeat_fraction_homopolymer() {
        // "AAAAAA" is 6 bases of poly-A run, "CGT" is 3 bases non-repeat → 6/9.
        let f = compute_repeat_fraction(b"AAAAAACGT");
        assert!(
            f > 0.5,
            "poly-A region should give repeat fraction > 0.5, got {f}"
        );
    }

    #[test]
    fn repeat_fraction_dinucleotide() {
        // ATATAT is 6 bases dinucleotide repeat.
        let f = compute_repeat_fraction(b"ATATATATCG");
        assert!(
            f > 0.5,
            "dinucleotide repeat should give fraction > 0.5, got {f}"
        );
    }

    #[test]
    fn repeat_fraction_no_repeats() {
        // ACGT has no run >= 3.
        let f = compute_repeat_fraction(b"ACGT");
        assert_eq!(f, 0.0, "ACGT has no repeats");
    }

    #[test]
    fn repeat_fraction_empty_is_zero() {
        assert_eq!(compute_repeat_fraction(b""), 0.0);
    }

    // ── Additional edge-case tests ───────────────────────────────────────────

    // Test: DUST score for a pure homopolymer run (64 A's) must be very high (> 25).
    // A 64-bp poly-A sequence has only one trinucleotide type (AAA), producing
    // the maximum possible score for that window size. The DUST formula yields
    // C(62,2) / 62 = (62 * 61 / 2) / 62 = 30.5 for a 64-bp homopolymer.
    #[test]
    fn dust_score_pure_homopolymer_very_high() {
        let seq = b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
        let score = compute_dust_score(seq, 64);
        assert!(
            score > 25.0,
            "pure 64-bp homopolymer must have DUST > 25, got {score}"
        );
    }

    // Test: DUST score for a dinucleotide repeat (ATATATATAT...) should be high.
    // Dinucleotide repeats are low-complexity and must score well above the
    // "borderline" threshold of 10.
    #[test]
    fn dust_score_dinucleotide_repeat_is_high() {
        let seq: Vec<u8> = b"AT".repeat(32).to_vec(); // 64 bp
        let score = compute_dust_score(&seq, 64);
        assert!(
            score > 10.0,
            "dinucleotide repeat (AT)x32 must have DUST > 10, got {score}"
        );
    }

    // Test: DUST score for a trinucleotide repeat (ACGACGACG...) should be moderate.
    // Trinucleotide repeats have 3 distinct trinucleotide types (ACG, CGA, GAC),
    // so the score is above random but below a pure homopolymer.
    #[test]
    fn dust_score_trinucleotide_repeat_moderate() {
        let seq: Vec<u8> = b"ACG".repeat(21).to_vec(); // 63 bp
        let score = compute_dust_score(&seq, 64);
        assert!(
            score > 3.0,
            "trinucleotide repeat must have noticeable DUST score, got {score}"
        );
        // Should be lower than a homopolymer.
        let homo_score = compute_dust_score(&[b'A'; 63], 63);
        assert!(
            score < homo_score,
            "trinucleotide DUST ({score}) must be less than homopolymer DUST ({homo_score})"
        );
    }

    // Test: DUST score for a very short sequence (< 3 bases) returns 0.0.
    // The algorithm requires trinucleotides; shorter inputs produce no windows.
    #[test]
    fn dust_score_short_sequence_returns_zero() {
        assert_eq!(
            compute_dust_score(b"AC", 64),
            0.0,
            "2-base sequence must return 0.0"
        );
        assert_eq!(
            compute_dust_score(b"A", 64),
            0.0,
            "1-base sequence must return 0.0"
        );
    }

    // Test: DUST score with N bases must not panic and must return a finite value.
    // Ambiguous bases (N) map to A (index 0) via base_index; the function must
    // remain numerically stable.
    #[test]
    fn dust_score_with_n_bases_does_not_panic() {
        let seq = b"NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN";
        let score = compute_dust_score(seq, 64);
        assert!(
            score.is_finite(),
            "DUST score with N bases must be finite, got {score}"
        );
    }

    // Test: repeat_fraction for all-same-base (AAAAAAA) returns 1.0.
    // A pure homopolymer of length >= 3 is entirely a repeat.
    #[test]
    fn repeat_fraction_all_same_base_is_one() {
        let f = compute_repeat_fraction(b"AAAAAAA");
        assert!(
            (f - 1.0).abs() < 1e-6,
            "pure homopolymer must have repeat fraction 1.0, got {f}"
        );
    }

    // Test: repeat_fraction for alternating dinucleotide repeat of length 8.
    // ACACACAC is a dinucleotide repeat of 4 units (8 bases, >= 6 threshold).
    #[test]
    fn repeat_fraction_dinucleotide_repeat_of_eight() {
        let f = compute_repeat_fraction(b"ACACACAC");
        assert!(
            f > 0.5,
            "AC dinucleotide repeat (8 bp) must have repeat fraction > 0.5, got {f}"
        );
    }

    // Test: repeat_fraction for a completely non-repetitive sequence.
    // ACGTTGCA has no homopolymer run >= 3 and no dinucleotide repeat >= 6.
    #[test]
    fn repeat_fraction_random_sequence_is_zero() {
        let f = compute_repeat_fraction(b"ACGTTGCA");
        assert_eq!(
            f, 0.0,
            "non-repetitive sequence must have repeat fraction 0.0, got {f}"
        );
    }

    // Test: repeat_fraction for a mixed sequence (AAAAACGT).
    // 5-base A run (>= 3 threshold) covers 5 of 8 bases = 0.625.
    #[test]
    fn repeat_fraction_mixed_sequence() {
        let f = compute_repeat_fraction(b"AAAAACGT");
        assert!(
            (f - 5.0 / 8.0).abs() < 1e-6,
            "AAAAACGT must have repeat fraction 5/8 = 0.625, got {f}"
        );
    }
}
