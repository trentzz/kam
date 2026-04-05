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

use crate::caller::{VariantCall, VariantType};

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
        let meta: ModelMeta = serde_json::from_slice(&meta_bytes)?;

        #[cfg(feature = "ml")]
        {
            let session = Session::builder()?
                .with_optimization_level(GraphOptimizationLevel::Level3)?
                .with_intra_threads(1)?
                .commit_from_file(model_path)?;
            Ok(Self { session, meta })
        }

        #[cfg(not(feature = "ml"))]
        {
            let _ = model_path;
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
    fn extract_features(&self, call: &VariantCall) -> Vec<f32> {
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

        let ref_alt_len_ratio = ref_len / (alt_len + 1.0);
        let indel_size = (ref_len - alt_len).abs();

        let duplex_enrichment = ndupalt / (vaf * alt_depth + 1e-9);
        let simplex_only_frac = nsimalt / (nalt + 1e-9);

        let conf_above_99 = if conf > 0.99 { 1.0_f32 } else { 0.0_f32 };
        let conf_above_999 = if conf > 0.999 { 1.0_f32 } else { 0.0_f32 };
        let sbp_above_05 = if sbp > 0.05 { 1.0_f32 } else { 0.0_f32 };

        let variant_class_enc = *self
            .meta
            .variant_class_map
            .get(variant_type_str(call.variant_type))
            .unwrap_or(&0) as f32;

        let lookup: HashMap<&str, f32> = [
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
