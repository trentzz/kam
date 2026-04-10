#!/usr/bin/env python3
"""Train and evaluate ML models on the twist-duplex dataset.

Loads training_data.csv produced by build_training_data_twist_duplex.py,
trains LightGBM and XGBoost on the rust-safe 33-feature set, and evaluates
on the held-out test split.

Outputs:
  docs/project/experiments/01-ml-twist-duplex/results/cv_results_twist_duplex.csv
  docs/project/experiments/01-ml-twist-duplex/results/test_results_twist_duplex.csv
  docs/project/experiments/01-ml-twist-duplex/results/feature_importance_twist_duplex.csv
  bigdata/experiments/01-ml-twist-duplex/models/lightgbm_twist_duplex.txt
  bigdata/experiments/01-ml-twist-duplex/models/xgboost_twist_duplex.json

Usage: python3 scripts/ml/train_eval_twist_duplex.py
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

warnings.filterwarnings("ignore", category=UndefinedMetricWarning)

import lightgbm as lgb
import xgboost as xgb

REPO = Path(__file__).resolve().parent.parent.parent
RESULTS = REPO / "docs/project/experiments/01-ml-twist-duplex/results"
MODELS  = REPO / "bigdata/experiments/01-ml-twist-duplex/models"
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

# Column aliases: the build script already renames these, but handle raw CSV
# column names defensively in case the CSV comes from an older run.
COL_ALIASES = {
    "n_molecules_ref": "nref",
    "n_molecules_alt": "nalt",
    "n_duplex_alt":    "ndupalt",
    "n_simplex_alt":   "nsimalt",
    "strand_bias_p":   "sbp",
    "confidence":      "conf",
    "vaf_ci_low":      "vaf_lo",
    "vaf_ci_high":     "vaf_hi",
}


def load_data(csv_path: Path) -> pd.DataFrame:
    """Load the training CSV, reading only columns needed for model training.

    Avoids loading large string columns (ref, alt, variant_type) which inflate
    pandas memory usage far beyond the on-disk file size.
    """
    needed = set(RUST_FEATURES) | {"label", "split", "sample_id"}
    for src in COL_ALIASES:
        needed.add(src)
    # Discover which needed columns actually exist in the file.
    header_cols = pd.read_csv(csv_path, nrows=0).columns.tolist()
    usecols = [c for c in header_cols if c in needed or COL_ALIASES.get(c) in needed]
    df = pd.read_csv(
        csv_path,
        usecols=usecols,
        dtype={c: "float32" for c in RUST_FEATURES if c in usecols},
    )
    df = df.rename(columns=COL_ALIASES)
    return df


def add_derived_features(df: pd.DataFrame) -> pd.DataFrame:
    """Add any derived features not already present in the CSV."""
    eps = 1e-9
    if "duplex_frac" not in df.columns:
        df["duplex_frac"] = df["ndupalt"] / (df["nalt"] + eps)
    if "has_duplex" not in df.columns:
        df["has_duplex"] = (df["ndupalt"] > 0).astype(float)
    if "ci_width" not in df.columns:
        lo_col = "vaf_lo" if "vaf_lo" in df.columns else "vaf_ci_low"
        hi_col = "vaf_hi" if "vaf_hi" in df.columns else "vaf_ci_high"
        df["ci_width"] = df[hi_col] - df[lo_col]
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
        df["ref_len"] = df["ref"].str.len() if "ref" in df.columns else 1.0
    if "alt_len" not in df.columns:
        df["alt_len"] = df["alt"].str.len() if "alt" in df.columns else 1.0
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
        "auroc":     roc_auc_score(y_true, y_pred_prob) if y_true.sum() > 0 else 0.0,
        "auprc":     average_precision_score(y_true, y_pred_prob) if y_true.sum() > 0 else 0.0,
        "precision": precision_score(y_true, y_pred, zero_division=0),
        "recall":    recall_score(y_true, y_pred, zero_division=0),
        "f1":        f1_score(y_true, y_pred, zero_division=0),
        "n_pos":     int(y_true.sum()),
        "n_total":   len(y_true),
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
    # Use the slimmed feature CSV (no ref/alt haplotype strings) if available.
    data_path = REPO / "bigdata/experiments/01-ml-twist-duplex/training_features.csv"
    if not data_path.exists():
        data_path = REPO / "bigdata/experiments/01-ml-twist-duplex/training_data.csv"
    if not data_path.exists():
        print(
            f"ERROR: training data not found.\n"
            "Run build_training_data_twist_duplex.py to generate it."
        )
        return
    print(f"Loading {data_path}...")
    df = load_data(data_path)
    df = add_derived_features(df)

    if "split" in df.columns:
        train_df = df[df["split"] == "train"].copy()
        test_df  = df[df["split"] == "test"].copy()
        print(f"Split found: {len(train_df)} train rows, {len(test_df)} test rows.")
    else:
        print("No 'split' column found; using all data as train (no held-out test).")
        train_df = df.copy()
        test_df  = pd.DataFrame()

    # Subsample negatives in the training set to keep memory and training time
    # manageable. Keep all positives + at most 10x negatives.
    neg_limit = int(train_df["label"].sum() * 10)
    neg_df = train_df[train_df["label"] == 0]
    pos_df = train_df[train_df["label"] == 1]
    if len(neg_df) > neg_limit:
        print(
            f"Subsampling negatives: {len(neg_df):,} → {neg_limit:,} "
            f"(keeping all {len(pos_df):,} positives)"
        )
        neg_df = neg_df.sample(neg_limit, random_state=42)
        train_df = pd.concat([pos_df, neg_df], ignore_index=True)
    else:
        print(f"No subsampling needed: {len(neg_df):,} negatives, {len(pos_df):,} positives.")

    available = [f for f in RUST_FEATURES if f in train_df.columns]
    missing   = [f for f in RUST_FEATURES if f not in train_df.columns]
    if missing:
        print(f"WARNING: {len(missing)} features missing from training data: {missing}")
    print(f"Using {len(available)} features: {available}")

    X_train = train_df[available].fillna(0.0)
    y_train = train_df["label"]
    groups  = train_df.get("sample_id", pd.Series(range(len(train_df))))

    pos_rate   = y_train.mean()
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
    lgb_rows = cross_validate(X_train, y_train, groups, "lightgbm_twist_duplex", lgb_clf)

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
    xgb_rows = cross_validate(X_train, y_train, groups, "xgboost_twist_duplex", xgb_clf)

    all_cv_rows = lgb_rows + xgb_rows
    cv_df  = pd.DataFrame(all_cv_rows)
    cv_path = RESULTS / "cv_results_twist_duplex.csv"
    cv_df.to_csv(cv_path, index=False)
    print(f"\nCV results → {cv_path}")

    # --- Final models on full training set ---
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

    lgb_model_path = MODELS / "lightgbm_twist_duplex.txt"
    xgb_model_path = MODELS / "xgboost_twist_duplex.json"
    lgb_final.booster_.save_model(str(lgb_model_path))
    xgb_final.save_model(str(xgb_model_path))
    print(f"Models saved: {lgb_model_path}, {xgb_model_path}")

    fi = pd.DataFrame({
        "feature":  available,
        "lgb_gain": lgb_final.booster_.feature_importance(importance_type="gain"),
        "xgb_gain": xgb_final.feature_importances_,
    }).sort_values("lgb_gain", ascending=False)
    fi_path = RESULTS / "feature_importance_twist_duplex.csv"
    fi.to_csv(fi_path, index=False)
    print(f"Feature importance → {fi_path}")

    # --- Held-out test evaluation ---
    if len(test_df) > 0:
        print("\nEvaluating on held-out test split...")
        X_test = test_df[available].fillna(0.0)
        y_test = test_df["label"]

        test_rows = []
        for name, clf in [("lightgbm_twist_duplex", lgb_final), ("xgboost_twist_duplex", xgb_final)]:
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
        test_path = RESULTS / "test_results_twist_duplex.csv"
        test_df_out.to_csv(test_path, index=False)
        print(f"Test results → {test_path}")
    else:
        print("\nNo held-out test data; skipping test evaluation.")

    # --- CV summary ---
    print("\nCV summary (mean ± std across 5 folds):")
    for model in ["lightgbm_twist_duplex", "xgboost_twist_duplex"]:
        sub = cv_df[cv_df["model"] == model]
        print(
            f"  {model}: AUPRC={sub['auprc'].mean():.4f}±{sub['auprc'].std():.4f}"
            f"  AUROC={sub['auroc'].mean():.4f}±{sub['auroc'].std():.4f}"
        )


if __name__ == "__main__":
    main()
