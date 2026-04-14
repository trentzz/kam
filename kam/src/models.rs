//! Built-in ML model registry.
//!
//! Resolves a model name (e.g. `"single-strand-v1"`) or a file path to an
//! [`MlScorer`]. Named models are bundled into the binary via `include_bytes!`.

use kam_call::ml::MlScorer;

static SINGLE_STRAND_V1_ONNX: &[u8] = include_bytes!("../models/single-strand-v1.onnx");
static SINGLE_STRAND_V1_META: &[u8] = include_bytes!("../models/single-strand-v1.json");

static TWIST_DUPLEX_V1_ONNX: &[u8] = include_bytes!("../models/twist-duplex-v1.onnx");
static TWIST_DUPLEX_V1_META: &[u8] = include_bytes!("../models/twist-duplex-v1.json");

static TWIST_DUPLEX_V2_ONNX: &[u8] = include_bytes!("../models/twist-duplex-v2.onnx");
static TWIST_DUPLEX_V2_META: &[u8] = include_bytes!("../models/twist-duplex-v2.json");

/// All built-in model names.
pub fn builtin_names() -> &'static [&'static str] {
    &["single-strand-v1", "twist-duplex-v1", "twist-duplex-v2"]
}

/// Resolve a model name or file path to an [`MlScorer`].
///
/// If `name_or_path` matches a built-in model name, load from bundled bytes.
/// Otherwise treat it as a file path: the companion metadata file must exist
/// at the same path with a `.json` extension replacing the original extension.
///
/// # Errors
///
/// Returns an error if the name is unknown or the files cannot be loaded.
pub fn resolve(name_or_path: &str) -> Result<MlScorer, Box<dyn std::error::Error>> {
    match name_or_path {
        "single-strand-v1" => MlScorer::from_bytes(SINGLE_STRAND_V1_ONNX, SINGLE_STRAND_V1_META),
        "twist-duplex-v1" => MlScorer::from_bytes(TWIST_DUPLEX_V1_ONNX, TWIST_DUPLEX_V1_META),
        "twist-duplex-v2" => MlScorer::from_bytes(TWIST_DUPLEX_V2_ONNX, TWIST_DUPLEX_V2_META),
        other => {
            let model_path = std::path::Path::new(other);
            let meta_path = model_path.with_extension("json");
            MlScorer::load(model_path, &meta_path)
        }
    }
}
