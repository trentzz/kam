"""Export a LightGBM classifier to ONNX for use with the Rust ML scorer.

Trains on the rust-safe 33-feature set derived from VariantCall fields,
then exports the model and companion metadata to bigdata/experiments/02-ml-single-strand/models/.

Usage:
    python scripts/ml/export_model.py

Requires: lightgbm, onnxmltools, skl2onnx, pandas, scikit-learn
"""

from __future__ import annotations

import json
import os
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
TRAINING_DATA = REPO_ROOT / "bigdata" / "experiments" / "02-ml-single-strand" / "training_data_v2.csv"
MODEL_DIR = REPO_ROOT / "bigdata" / "experiments" / "02-ml-single-strand" / "models"
MODEL_PATH = MODEL_DIR / "lightgbm_rust.onnx"
META_PATH = MODEL_DIR / "model_meta.json"

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
    if not TRAINING_DATA.exists():
        print(
            f"Training data not found at {TRAINING_DATA}. "
            "No model exported — rerun after generating training data."
        )
        sys.exit(0)

    try:
        import lightgbm as lgb
        import numpy as np
        import pandas as pd
        from onnxmltools import convert_lightgbm
        from onnxmltools.convert.common.data_types import FloatTensorType
    except ImportError as exc:
        print(f"Missing dependency: {exc}. Install lightgbm, onnxmltools, and skl2onnx.")
        sys.exit(1)

    print(f"Loading training data from {TRAINING_DATA} ...")
    df = pd.read_csv(TRAINING_DATA)

    # Keep only rows that have all required feature columns.
    available = [c for c in RUST_FEATURES if c in df.columns]
    missing = [c for c in RUST_FEATURES if c not in df.columns]
    if missing:
        print(f"WARNING: {len(missing)} feature(s) not in training data: {missing}")
        print("Proceeding with available features only.")

    if "label" not in df.columns:
        print("ERROR: 'label' column not found in training data.")
        sys.exit(1)

    df = df.dropna(subset=available + ["label"])
    X = df[available].values.astype(np.float32)
    y = df["label"].values.astype(int)

    print(f"Training on {len(X)} samples, {len(available)} features ...")

    clf = lgb.LGBMClassifier(
        is_unbalance=True,
        n_estimators=200,
        max_depth=6,
        num_leaves=31,
        random_state=42,
    )
    clf.fit(X, y)

    print("Exporting to ONNX ...")
    MODEL_DIR.mkdir(parents=True, exist_ok=True)

    initial_types = [("X", FloatTensorType([None, len(available)]))]
    onnx_model = convert_lightgbm(clf.booster_, initial_types=initial_types)

    with open(MODEL_PATH, "wb") as f:
        f.write(onnx_model.SerializeToString())
    print(f"ONNX model saved to {MODEL_PATH}")

    meta = {
        "version": "1",
        "feature_names": available,
        "ml_pass_threshold": ML_PASS_THRESHOLD,
        "variant_class_map": VARIANT_CLASS_MAP,
    }
    with open(META_PATH, "w") as f:
        json.dump(meta, f, indent=2)
    print(f"Model metadata saved to {META_PATH}")


if __name__ == "__main__":
    main()
