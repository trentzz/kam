"""Export the trained twist-duplex LightGBM model to ONNX.

Loads lightgbm_twist_duplex.txt produced by train_eval_twist_duplex.py
and exports it to ONNX for use with the Rust ML scorer. Does not retrain.

Usage:
    python scripts/ml/export_model_twist_duplex.py

Requires: lightgbm, onnxmltools, skl2onnx
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

REPO_ROOT    = Path(__file__).resolve().parents[2]
MODEL_DIR    = REPO_ROOT / "bigdata/experiments/01-ml-twist-duplex/models"
SOURCE_MODEL = MODEL_DIR / "lightgbm_twist_duplex.txt"
MODEL_PATH   = MODEL_DIR / "lightgbm_rust.onnx"
META_PATH    = MODEL_DIR / "model_meta.json"

# Features computable from VariantCall at inference time (no pipeline state needed).
RUST_FEATURES = [
    "vaf", "nref", "nalt", "ndupalt", "nsimalt", "sbp", "conf",
    "ref_len", "alt_len",
    "duplex_frac", "has_duplex", "ci_width", "alt_depth",
    "log_nalt", "log_nref", "log_alt_depth", "log_vaf",
    "vaf_times_conf", "vaf_times_nalt",
    "nalt_over_conf", "ci_width_rel", "snr",
    "conf_sq", "nalt_sq", "vaf_sq",
    "ref_alt_len_ratio", "indel_size",
    "duplex_enrichment", "simplex_only_frac",
    "conf_above_99", "conf_above_999", "sbp_above_05",
    "variant_class_enc",
]

VARIANT_CLASS_MAP = {
    "SNV": 0,
    "Insertion": 1,
    "Deletion": 2,
    "MNV": 3,
    "Complex": 4,
    "LargeDeletion": 5,
    "TandemDuplication": 6,
    "Inversion": 7,
    "Fusion": 8,
    "InvDel": 9,
    "NovelInsertion": 10,
}

ML_PASS_THRESHOLD = 0.5


def main() -> None:
    if not SOURCE_MODEL.exists():
        print(
            f"Trained model not found at {SOURCE_MODEL}. "
            "Run train_eval_twist_duplex.py first to produce it."
        )
        sys.exit(1)

    try:
        import lightgbm as lgb
        from onnxmltools import convert_lightgbm
        from onnxmltools.convert.common.data_types import FloatTensorType
    except ImportError as exc:
        print(f"Missing dependency: {exc}. Install lightgbm, onnxmltools, and skl2onnx.")
        sys.exit(1)

    print(f"Loading trained model from {SOURCE_MODEL} ...")
    booster = lgb.Booster(model_file=str(SOURCE_MODEL))

    n_features = len(RUST_FEATURES)
    print(f"Exporting {n_features}-feature model to ONNX ...")
    MODEL_DIR.mkdir(parents=True, exist_ok=True)

    initial_types = [("X", FloatTensorType([None, n_features]))]
    onnx_model = convert_lightgbm(booster, initial_types=initial_types)

    with open(MODEL_PATH, "wb") as f:
        f.write(onnx_model.SerializeToString())
    print(f"ONNX model saved to {MODEL_PATH}")

    meta = {
        "version":           "twist-duplex-1",
        "feature_names":     RUST_FEATURES,
        "ml_pass_threshold": ML_PASS_THRESHOLD,
        "variant_class_map": VARIANT_CLASS_MAP,
    }
    with open(META_PATH, "w") as f:
        json.dump(meta, f, indent=2)
    print(f"Model metadata saved to {META_PATH}")


if __name__ == "__main__":
    main()
