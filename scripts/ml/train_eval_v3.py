#!/usr/bin/env python3
"""Train and evaluate ML models on the ML3 dataset with an explicit train/test split.

This script uses training_data_v3.csv, which combines data from all three
generations of varforge configs (original 250 + ML1 433 + ML2 3433 + ML3 train 10000).
The ML3 test split (1000 configs) is held out entirely until final evaluation.

Workflow:
  1. Load training data, filter to the ML3 train split.
  2. Train LightGBM and XGBoost on the rust-safe 33-feature set.
  3. Evaluate on the held-out ML3 test split.
  4. Export the best model to ONNX via export_model.py.

Outputs:
  docs/benchmarking/ml/results/cv_results_v3.csv         — per-fold CV metrics
  docs/benchmarking/ml/results/test_results_v3.csv       — held-out test metrics
  docs/benchmarking/ml/results/feature_importance_v3.csv — feature gain
  docs/benchmarking/ml/models/lightgbm_v3.txt            — trained LightGBM
  docs/benchmarking/ml/models/xgboost_v3.json            — trained XGBoost

Usage: python3 scripts/ml/train_eval_v3.py
Run from the repository root.
"""

from __future__ import annotations

import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.exceptions import UndefinedMetricWarning
from sklearn.metrics import (
    average_precision_score,
    precision_score,
    recall_score,
    f1_score,
    roc_auc_score,
)
from sklearn.model_selection import GroupKFold
from sklearn.preprocessing import LabelEncoder

warnings.filterwarnings("ignore", category=UndefinedMetricWarning)

import lightgbm as lgb
import xgboost as xgb

REPO = Path(__file__).resolve().parent.parent.parent
RESULTS = REPO / "docs/benchmarking/ml/results"
MODELS = REPO / "docs/benchmarking/ml/models"
RESULTS.mkdir(parents=True, exist_ok=True)
MODELS.mkdir(parents=True, exist_ok=True)

# Rust-safe features: all directly computable from VariantCall at inference time.
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

# Column aliases from the training data CSV to the RUST_FEATURES names.
COL_ALIASES = {
    "n_molecules_ref": "nref",
    "n_molecules_alt": "nalt",
    "n_duplex_alt": "ndupalt",
    "n_simplex_alt": "nsimalt",
    "strand_bias_p": "sbp",
    "confidence": "conf",
}


def load_data(csv_path: Path) -> pd.DataFrame:
    """Load and preprocess the training CSV."""
    df = pd.read_csv(csv_path)
    df = df.rename(columns=COL_ALIASES)
    return df


def add_derived_features(df: pd.DataFrame) -> pd.DataFrame:
    """Add any derived features not already in the CSV."""
    eps = 1e-9
    if "duplex_frac" not in df.columns:
        df["duplex_frac"] = df["ndupalt"] / (df["nalt"] + eps)
    if "has_duplex" not in df.columns:
        df["has_duplex"] = (df["ndupalt"] > 0).astype(float)
    if "ci_width" not in df.columns:
        df["ci_width"] = df["vaf_ci_high"] - df["vaf_ci_low"]
    if "alt_depth" not in df.columns:
        df["alt_depth"] = df["nref"] + df["nalt"]
    if "log_nalt" not in df.columns:
        df["log_nalt"] = np.log(df["nalt"] + 1)
    if "log_nref" not in df.columns:
        df["log_nref"] = np.log(df["nref"] + 1)
    if "log_alt_depth" not in df.columns:
        df["log_alt_depth"] = np.log(df["alt_depth"] + 1)
    if "log_vaf" not in df.columns:
        df["log_vaf"] = np.log(df["vaf"] + 1e-6)
    if "vaf_times_conf" not in df.columns:
        df["vaf_times_conf"] = df["vaf"] * df["conf"]
    if "vaf_times_nalt" not in df.columns:
        df["vaf_times_nalt"] = df["vaf"] * df["nalt"]
    if "nalt_over_conf" not in df.columns:
        df["nalt_over_conf"] = df["nalt"] / (df["conf"] + eps)
    if "ci_width_rel" not in df.columns:
        df["ci_width_rel"] = df["ci_width"] / (df["vaf"] + eps)
    if "snr" not in df.columns:
        df["snr"] = df["nalt"] / (df["nref"] + 1)
    if "conf_sq" not in df.columns:
        df["conf_sq"] = df["conf"] ** 2
    if "nalt_sq" not in df.columns:
        df["nalt_sq"] = df["nalt"] ** 2
    if "vaf_sq" not in df.columns:
        df["vaf_sq"] = df["vaf"] ** 2
    if "ref_len" not in df.columns:
        df["ref_len"] = df["ref_seq"].str.len() if "ref_seq" in df.columns else 1.0
    if "alt_len" not in df.columns:
        df["alt_len"] = df["alt_seq"].str.len() if "alt_seq" in df.columns else 1.0
    if "ref_alt_len_ratio" not in df.columns:
        df["ref_alt_len_ratio"] = df["ref_len"] / (df["alt_len"] + 1)
    if "indel_size" not in df.columns:
        df["indel_size"] = (df["ref_len"] - df["alt_len"]).abs()
    if "duplex_enrichment" not in df.columns:
        df["duplex_enrichment"] = df["ndupalt"] / (df["vaf"] * df["alt_depth"] + eps)
    if "simplex_only_frac" not in df.columns:
        df["simplex_only_frac"] = df["nsimalt"] / (df["nalt"] + eps)
    if "conf_above_99" not in df.columns:
        df["conf_above_99"] = (df["conf"] > 0.99).astype(float)
    if "conf_above_999" not in df.columns:
        df["conf_above_999"] = (df["conf"] > 0.999).astype(float)
    if "sbp_above_05" not in df.columns:
        df["sbp_above_05"] = (df["sbp"] > 0.05).astype(float)

    # Encode variant type if not yet done.
    if "variant_class_enc" not in df.columns and "variant_type" in df.columns:
        vt_map = {
            "SNV": 0, "Insertion": 1, "Deletion": 2, "MNV": 3,
            "Complex": 4, "LargeDeletion": 5, "TandemDuplication": 6,
            "Inversion": 7, "Fusion": 8, "InvDel": 9, "NovelInsertion": 10,
        }
        df["variant_class_enc"] = df["variant_type"].map(vt_map).fillna(0)

    return df


def evaluate(y_true, y_pred_prob, threshold: float = 0.5) -> dict:
    """Compute classification metrics at a given threshold."""
    y_pred = (y_pred_prob >= threshold).astype(int)
    return {
        "auroc": roc_auc_score(y_true, y_pred_prob) if y_true.sum() > 0 else 0.0,
        "auprc": average_precision_score(y_true, y_pred_prob) if y_true.sum() > 0 else 0.0,
        "precision": precision_score(y_true, y_pred, zero_division=0),
        "recall": recall_score(y_true, y_pred, zero_division=0),
        "f1": f1_score(y_true, y_pred, zero_division=0),
        "n_pos": int(y_true.sum()),
        "n_total": len(y_true),
    }


def cross_validate(X, y, groups, model_name: str, clf) -> list[dict]:
    """5-fold GroupKFold cross-validation. Returns per-fold metric dicts."""
    gkf = GroupKFold(n_splits=5)
    rows = []
    for fold, (train_idx, val_idx) in enumerate(gkf.split(X, y, groups)):
        X_tr, X_va = X.iloc[train_idx], X.iloc[val_idx]
        y_tr, y_va = y.iloc[train_idx], y.iloc[val_idx]
        clf.fit(X_tr, y_tr)
        probs = clf.predict_proba(X_va)[:, 1]
        m = evaluate(y_va, probs)
        m.update({"model": model_name, "fold": fold})
        rows.append(m)
        print(
            f"  fold {fold}: AUPRC={m['auprc']:.4f}  AUROC={m['auroc']:.4f}"
            f"  P={m['precision']:.4f}  R={m['recall']:.4f}"
        )
    return rows


def main() -> None:
    """Main training and evaluation loop."""
    # Try loading v3 data first; fall back to v2.
    data_v3 = REPO / "docs/benchmarking/ml/training_data_v3.csv"
    data_v2 = REPO / "docs/benchmarking/ml/training_data_v2.csv"
    if data_v3.exists():
        print(f"Loading {data_v3}...")
        df = load_data(data_v3)
    elif data_v2.exists():
        print(f"training_data_v3.csv not found; loading {data_v2}...")
        df = load_data(data_v2)
    else:
        print("No training data found. Run build_training_data_v2.py first.")
        return

    df = add_derived_features(df)

    # Identify split column: 'split' (ml3) or fall back to all as train.
    if "split" in df.columns:
        train_df = df[df["split"] == "train"].copy()
        test_df = df[df["split"] == "test"].copy()
        print(
            f"Split found: {len(train_df)} train rows, {len(test_df)} test rows."
        )
    else:
        print("No 'split' column found; using all data as train (no held-out test).")
        train_df = df.copy()
        test_df = pd.DataFrame()

    # Check available features.
    available = [f for f in RUST_FEATURES if f in train_df.columns]
    missing = [f for f in RUST_FEATURES if f not in train_df.columns]
    if missing:
        print(f"WARNING: {len(missing)} features missing from training data: {missing}")
    print(f"Using {len(available)} features: {available}")

    X_train = train_df[available].fillna(0.0)
    y_train = train_df["label"]
    groups = train_df.get("sample_id", pd.Series(range(len(train_df))))

    pos_rate = y_train.mean()
    pos_weight = (1 - pos_rate) / (pos_rate + 1e-9)
    print(
        f"Training set: {len(y_train)} rows, {int(y_train.sum())} positives "
        f"({pos_rate:.3%}), scale_pos_weight={pos_weight:.1f}"
    )

    # --- LightGBM ---
    print("\nLightGBM 5-fold CV...")
    lgb_clf = lgb.LGBMClassifier(
        n_estimators=300,
        max_depth=6,
        num_leaves=31,
        learning_rate=0.05,
        is_unbalance=True,
        random_state=42,
        n_jobs=4,
        verbosity=-1,
    )
    lgb_rows = cross_validate(X_train, y_train, groups, "lightgbm_v3", lgb_clf)

    # --- XGBoost ---
    print("\nXGBoost 5-fold CV...")
    xgb_clf = xgb.XGBClassifier(
        n_estimators=300,
        max_depth=6,
        learning_rate=0.05,
        scale_pos_weight=pos_weight,
        eval_metric="aucpr",
        random_state=42,
        n_jobs=4,
        verbosity=0,
    )
    xgb_rows = cross_validate(X_train, y_train, groups, "xgboost_v3", xgb_clf)

    all_cv_rows = lgb_rows + xgb_rows
    cv_df = pd.DataFrame(all_cv_rows)
    cv_path = RESULTS / "cv_results_v3.csv"
    cv_df.to_csv(cv_path, index=False)
    print(f"\nCV results → {cv_path}")

    # --- Final model on full training set ---
    print("\nFitting final models on full training set...")
    lgb_final = lgb.LGBMClassifier(
        n_estimators=300,
        max_depth=6,
        num_leaves=31,
        learning_rate=0.05,
        is_unbalance=True,
        random_state=42,
        n_jobs=4,
        verbosity=-1,
    )
    lgb_final.fit(X_train, y_train)

    xgb_final = xgb.XGBClassifier(
        n_estimators=300,
        max_depth=6,
        learning_rate=0.05,
        scale_pos_weight=pos_weight,
        eval_metric="aucpr",
        random_state=42,
        n_jobs=4,
        verbosity=0,
    )
    xgb_final.fit(X_train, y_train)

    # Save models.
    lgb_model_path = MODELS / "lightgbm_v3.txt"
    xgb_model_path = MODELS / "xgboost_v3.json"
    lgb_final.booster_.save_model(str(lgb_model_path))
    xgb_final.save_model(str(xgb_model_path))
    print(f"Models saved: {lgb_model_path}, {xgb_model_path}")

    # Feature importance.
    fi = pd.DataFrame({
        "feature": available,
        "lgb_gain": lgb_final.booster_.feature_importance(importance_type="gain"),
        "xgb_gain": xgb_final.feature_importances_,
    }).sort_values("lgb_gain", ascending=False)
    fi_path = RESULTS / "feature_importance_v3.csv"
    fi.to_csv(fi_path, index=False)
    print(f"Feature importance → {fi_path}")

    # --- Held-out test evaluation ---
    if len(test_df) > 0:
        print("\nEvaluating on held-out test split...")
        X_test = test_df[available].fillna(0.0)
        y_test = test_df["label"]

        test_rows = []
        for name, clf in [("lightgbm_v3", lgb_final), ("xgboost_v3", xgb_final)]:
            probs = clf.predict_proba(X_test)[:, 1]
            m = evaluate(y_test, probs)
            m["model"] = name
            test_rows.append(m)
            print(
                f"  {name}: AUPRC={m['auprc']:.4f}  AUROC={m['auroc']:.4f}"
                f"  P={m['precision']:.4f}  R={m['recall']:.4f}"
                f"  n_pos={m['n_pos']}/{m['n_total']}"
            )

        test_df_out = pd.DataFrame(test_rows)
        test_path = RESULTS / "test_results_v3.csv"
        test_df_out.to_csv(test_path, index=False)
        print(f"Test results → {test_path}")
    else:
        print("\nNo held-out test data; skipping test evaluation.")

    # --- CV summary ---
    print("\nCV summary (mean ± std across 5 folds):")
    for model in ["lightgbm_v3", "xgboost_v3"]:
        sub = cv_df[cv_df["model"] == model]
        print(
            f"  {model}: AUPRC={sub['auprc'].mean():.4f}±{sub['auprc'].std():.4f}"
            f"  AUROC={sub['auroc'].mean():.4f}±{sub['auroc'].std():.4f}"
        )


if __name__ == "__main__":
    main()
