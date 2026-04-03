# ML Integration Design

## Context

`kam call` outputs a TSV with one row per variant call. This document covers how to add ML-based post-call filtering. The model receives per-call features and outputs a probability score and a filter label. The design question is how to ship the model and run inference inside the Rust binary.

---

## Options

### Option 1: `lightgbm` Rust crate

The `lightgbm` crate wraps the LightGBM C API via FFI. It loads a `.txt` model file at runtime and runs inference in process.

**Pros**

- No subprocess overhead.
- Single binary deployment; model file ships alongside.
- Full control over inference latency.
- LightGBM C inference is fast and deterministic.

**Cons**

- Requires linking against the LightGBM C library (`libLightGBM.so` or static). This complicates cross-compilation and CI.
- The `lightgbm` crate is not maintained as actively as the Python library. API may lag behind model format changes.
- Builds become platform-specific unless the static library is vendored.

---

### Option 2: ONNX export with `ort`

Export the trained LightGBM model to ONNX using `sklearn2pmml` or `onnxmltools`. Load it at runtime with the `ort` crate (ONNX Runtime Rust bindings).

**Pros**

- ONNX Runtime is well-maintained and supports many backends (CPU, GPU, etc.).
- Model format is framework-agnostic. Switching from LightGBM to another model type requires no Rust code change.
- The `ort` crate has good ergonomics and active maintenance.

**Cons**

- Adds a dependency on ONNX Runtime, which is a large C++ library (~50 MB static).
- LightGBM-to-ONNX conversion is not always lossless. Tree structure and feature handling must be verified.
- Two-step process: train in Python, export to ONNX, then test the exported model matches the original.

---

### Option 3: Python subprocess at call time

Invoke a Python script from within `kam call` using `std::process::Command`. Pass variant features via stdin (JSON or TSV) and read probabilities from stdout.

**Pros**

- No FFI or linking required.
- Model can be any Python-supported framework. Updating the model requires no Rust rebuild.

**Cons**

- Python must be available on the execution host. This is not guaranteed in HPC environments.
- Subprocess startup cost is 100-300 ms per invocation. At call time this is called per-sample, so the total overhead is one startup, but it still complicates the binary.
- Harder to make deterministic and hermetic. Environment differences (Python version, package versions) can change output.
- Piping TSV/JSON adds serialisation complexity.

---

## Recommendation: ONNX with `ort`

Use ONNX export and the `ort` crate.

The key reason is flexibility without coupling. ONNX Runtime handles inference correctly and efficiently without requiring a live Python environment. The `ort` crate supports loading a model file at runtime, so the model can be updated without recompiling. LightGBM models export to ONNX reliably for the tree-based inference path used here.

The binary size cost is real but acceptable. Rust builds already link large libraries (Rayon, DashMap, etc.). Linking ONNX Runtime statically adds roughly 50 MB but produces a self-contained binary that runs on any x86-64 Linux host, which is the target deployment environment.

The `lightgbm` Rust crate is second choice. Its FFI approach works but the maintenance risk is higher and cross-compilation is harder. The subprocess approach is unsuitable for production use.

---

## CLI Changes

Add one optional flag to `kam call`:

```
--ml-model <path>    Path to ONNX model file. If omitted, ML scoring is skipped.
```

When `--ml-model` is provided:

1. Load the model at startup using `ort::Session::builder().commit_from_file(path)`.
2. For each call, assemble the feature vector in the same order as training.
3. Run inference. Read the positive-class probability from the output tensor.
4. Write two additional columns to the TSV output.

---

## Output Changes

Add two columns to the variant TSV, appended after existing columns:

| Column | Type | Description |
|--------|------|-------------|
| `ml_prob` | f32 | Probability that the call is a true variant. Range [0, 1]. `NaN` if no model was loaded. |
| `ml_filter` | string | `PASS` if `ml_prob >= 0.5`, `MLFilt` otherwise. Empty string if no model loaded. |

These columns appear only when `--ml-model` is given. Downstream tools that do not provide the flag receive no change to the TSV schema.

---

## Model Loading

```rust
use ort::{Environment, Session, SessionBuilder};

pub struct MlScorer {
    session: Session,
}

impl MlScorer {
    pub fn load(path: &std::path::Path) -> anyhow::Result<Self> {
        let environment = Environment::builder().with_name("kam_ml").build()?.into_arc();
        let session = SessionBuilder::new(&environment)?
            .with_intra_threads(1)?
            .commit_from_file(path)?;
        Ok(Self { session })
    }

    /// Score one call. Returns probability in [0, 1].
    pub fn score(&self, features: &[f32]) -> anyhow::Result<f32> {
        use ort::inputs;
        let input = ndarray::Array2::from_shape_vec((1, features.len()), features.to_vec())?;
        let outputs = self.session.run(inputs!["float_input" => input.view()]?)?;
        // LightGBM ONNX output: [n, 2] probabilities, column 1 is positive class
        let probs = outputs[1].try_extract_tensor::<f32>()?;
        Ok(probs[[0, 1]])
    }
}
```

The feature order must match training exactly. Define a `const FEATURE_ORDER: &[&str]` in the Rust source that mirrors `FEATURE_COLS` in `train_eval.py`. Categorical features require the same integer encoding used during training; store the encoding maps alongside the ONNX model or embed them as a small JSON sidecar.

---

## Inference Integration

Inference runs per call in the output loop of `kam call`. The overhead is negligible for tree-based models (microseconds per call). No batching is required for the call counts typical in a panel sample (thousands to tens of thousands of calls).

---

## Model Update Workflow

1. Retrain in Python using `train_eval.py`.
2. Export best model to ONNX: `booster.save_model("model.txt")` then convert with `onnxmltools`.
3. Validate ONNX output matches Python output on a held-out sample.
4. Replace the ONNX file on the deployment host. No binary rebuild needed.
