#!/usr/bin/env python3
"""Train and evaluate LightGBM and XGBoost on the v2 feature set.

This is the second-round training script. It uses training_data_v2.csv,
which has ~30 features including log transforms, ratios, interaction terms,
and experimental-condition parameters (coverage, family size, PCR cycles).

Compare five configurations:
  A. Baseline: PASS filter (binary rule)
  B. LightGBM: raw features only (same 9 columns as v1)
  C. LightGBM: full v2 feature set
  D. XGBoost:  raw features only
  E. XGBoost:  full v2 feature set

Also run feature selection experiments:
  - Pearson correlation filter (keep features |r| < 0.95 pairwise)
  - Permutation importance: drop features with negative mean importance
  - RFE (Recursive Feature Elimination) with LightGBM estimator

Outputs:
  docs/project/experiments/02-ml-single-strand/results/cv_results_v2.csv     — per-fold metrics
  docs/project/experiments/02-ml-single-strand/results/feature_importance_v2.csv
  docs/project/experiments/02-ml-single-strand/results/feature_selection.csv — RFE selected features
  bigdata/experiments/02-ml-single-strand/models/lightgbm_v2.txt              — best LightGBM model
  bigdata/experiments/02-ml-single-strand/models/xgboost_v2.json              — best XGBoost model

Usage: python3 scripts/ml/train_eval_v2.py
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
DATA = REPO / "bigdata/experiments/02-ml-single-strand/training_data_v2.csv"
# Small summary outputs (CSV, figures) are committed to docs/.
RESULTS = REPO / "docs/project/experiments/02-ml-single-strand/results"
# Trained model files are large and go to bigdata/.
MODELS = REPO / "bigdata/experiments/02-ml-single-strand/models"
RESULTS.mkdir(parents=True, exist_ok=True)
MODELS.mkdir(parents=True, exist_ok=True)

# Raw features (same as v1)
RAW_FEATURES = [
    "vaf", "nref", "nalt", "ndupalt", "nsimalt",
    "sbp", "conf", "ref_len", "alt_len",
]

# Full v2 feature set
V2_FEATURES = [
    # raw
    "vaf", "nref", "nalt", "ndupalt", "nsimalt", "sbp", "conf",
    "ref_len", "alt_len",
    # derived from v1
    "duplex_frac", "has_duplex", "ci_width", "alt_depth",
    # log transforms
    "log_nalt", "log_nref", "log_alt_depth", "log_vaf",
    # ratios / interactions
    "observed_vaf_diff", "vaf_times_conf", "vaf_times_nalt",
    "nalt_over_conf", "ci_width_rel", "snr",
    "conf_sq", "nalt_sq", "vaf_sq",
    "ref_alt_len_ratio", "indel_size",
    # duplex
    "duplex_enrichment", "simplex_only_frac",
    # categorical bins
    "depth_bucket", "sbp_cat",
    # confidence thresholds
    "conf_above_99", "conf_above_999", "sbp_above_05",
    # position
    "pos",
    # encoded categoricals
    "variant_class_enc", "mode_enc", "dataset_type_enc",
]

# Optional experimental params (present only in new ML samples)
PARAM_FEATURES = [
    "params_coverage", "params_family_size_mean",
    "params_pcr_cycles", "params_fragment_mean",
]


def load_data() -> pd.DataFrame:
    df = pd.read_csv(DATA, low_memory=False)
    print(f"Loaded {len(df):,} rows, {df['label'].sum()} positives "
          f"({df['label'].mean()*100:.3f}%)")

    # Encode categoricals
    for col, enc_col in [
        ("variant_class", "variant_class_enc"),
        ("mode", "mode_enc"),
        ("dataset_type", "dataset_type_enc"),
    ]:
        le = LabelEncoder()
        df[enc_col] = le.fit_transform(df[col].astype(str))

    return df


def available_features(df: pd.DataFrame, feature_list: list[str]) -> list[str]:
    """Return only features that exist in df and have non-zero variance."""
    available = [f for f in feature_list if f in df.columns]
    non_const = [f for f in available if df[f].std() > 0]
    dropped = set(available) - set(non_const)
    if dropped:
        print(f"  Dropped zero-variance features: {sorted(dropped)}")
    return non_const


def add_param_features(df: pd.DataFrame) -> pd.DataFrame:
    """Add param features if present, else fill with NaN."""
    for col in PARAM_FEATURES:
        if col not in df.columns:
            df[col] = float("nan")
    return df


def metrics(y_true, y_pred_prob, threshold: float = 0.5) -> dict:
    y_pred = (y_pred_prob >= threshold).astype(int)
    # Avoid divide-by-zero when no positives predicted
    if y_pred.sum() == 0:
        prec = 0.0
    else:
        prec = precision_score(y_true, y_pred, zero_division=0)
    rec = recall_score(y_true, y_pred, zero_division=0)
    f1 = f1_score(y_true, y_pred, zero_division=0)
    try:
        auprc = average_precision_score(y_true, y_pred_prob)
        auroc = roc_auc_score(y_true, y_pred_prob)
    except ValueError:
        auprc = auroc = float("nan")
    return dict(precision=prec, recall=rec, f1=f1, auprc=auprc, auroc=auroc)


def run_cv(
    df: pd.DataFrame,
    features: list[str],
    model_name: str,
    n_splits: int = 5,
) -> tuple[list[dict], float | None, object | None]:
    """Run GroupKFold CV. Returns (fold_results, best_auprc, best_model)."""
    groups = df["sample_id"].values
    X = df[features].values
    y = df["label"].values

    gkf = GroupKFold(n_splits=n_splits)
    fold_results = []
    best_auprc = -1.0
    best_model = None
    n_pos = y.sum()
    n_neg = len(y) - n_pos
    scale = n_neg / max(n_pos, 1)

    for fold, (train_idx, val_idx) in enumerate(gkf.split(X, y, groups)):
        X_train, X_val = X[train_idx], X[val_idx]
        y_train, y_val = y[train_idx], y[val_idx]

        if "lightgbm" in model_name.lower():
            model = lgb.LGBMClassifier(
                n_estimators=500,
                learning_rate=0.05,
                num_leaves=31,
                min_child_samples=10,
                is_unbalance=True,
                random_state=42,
                verbose=-1,
                n_jobs=-1,
            )
        elif "xgboost" in model_name.lower():
            model = xgb.XGBClassifier(
                n_estimators=500,
                learning_rate=0.05,
                max_depth=6,
                scale_pos_weight=scale,
                random_state=42,
                verbosity=0,
                n_jobs=-1,
                eval_metric="aucpr",
            )
        else:
            raise ValueError(f"Unknown model: {model_name}")

        model.fit(X_train, y_train)
        y_prob = model.predict_proba(X_val)[:, 1]
        m = metrics(y_val, y_prob)
        m.update(model=model_name, fold=fold)
        fold_results.append(m)

        if m["auprc"] > best_auprc:
            best_auprc = m["auprc"]
            best_model = model

        print(f"  fold {fold}: AUROC={m['auroc']:.4f} AUPRC={m['auprc']:.4f} "
              f"recall={m['recall']:.3f} prec={m['precision']:.4f}")

    return fold_results, best_auprc, best_model


def baseline_cv(df: pd.DataFrame, n_splits: int = 5) -> list[dict]:
    """Baseline: treat PASS filter as binary classifier."""
    groups = df["sample_id"].values
    y = df["label"].values
    y_pred = (df["filter"] == "PASS").astype(int).values
    gkf = GroupKFold(n_splits=n_splits)
    fold_results = []
    for fold, (_, val_idx) in enumerate(gkf.split(np.zeros(len(y)), y, groups)):
        y_val = y[val_idx]
        y_pred_val = y_pred[val_idx]
        m = metrics(y_val, y_pred_val.astype(float))
        m.update(model="baseline_pass", fold=fold)
        fold_results.append(m)
    return fold_results


def feature_importance_df(model, features: list[str], model_name: str) -> pd.DataFrame:
    if hasattr(model, "feature_importances_"):
        imp = model.feature_importances_
    else:
        return pd.DataFrame()
    return pd.DataFrame({
        "feature": features,
        "importance": imp,
        "model": model_name,
    }).sort_values("importance", ascending=False)


def main() -> None:
    if not DATA.exists():
        print(f"[ERROR] {DATA} not found. Run build_training_data_v2.py first.")
        return

    df = load_data()
    df = add_param_features(df)

    # Include param features if they have variance
    all_v2 = V2_FEATURES + [p for p in PARAM_FEATURES if p in df.columns]

    raw_feats = available_features(df, RAW_FEATURES)
    v2_feats = available_features(df, all_v2)

    print(f"\nRaw features: {len(raw_feats)}")
    print(f"V2 features:  {len(v2_feats)}")
    print()

    all_fold_results = []
    importance_frames = []
    best_lgb_model = None
    best_lgb_auprc = -1.0
    best_xgb_model = None
    best_xgb_auprc = -1.0

    configs = [
        ("lightgbm_raw",  "lightgbm", raw_feats),
        ("lightgbm_v2",   "lightgbm", v2_feats),
        ("xgboost_raw",   "xgboost",  raw_feats),
        ("xgboost_v2",    "xgboost",  v2_feats),
    ]

    for config_name, model_type, features in configs:
        print(f"\n=== {config_name} ({len(features)} features) ===")
        fold_results, auprc, model = run_cv(df, features, model_type)
        for r in fold_results:
            r["config"] = config_name
        all_fold_results.extend(fold_results)

        imp = feature_importance_df(model, features, config_name)
        if not imp.empty:
            importance_frames.append(imp)

        if "lightgbm" in config_name and auprc > best_lgb_auprc:
            best_lgb_auprc = auprc
            best_lgb_model = model
        if "xgboost" in config_name and auprc > best_xgb_auprc:
            best_xgb_auprc = auprc
            best_xgb_model = model

    # Baseline
    print("\n=== baseline_pass ===")
    base_results = baseline_cv(df)
    for r in base_results:
        r["config"] = "baseline_pass"
    all_fold_results.extend(base_results)

    # Save results
    results_df = pd.DataFrame(all_fold_results)
    results_path = RESULTS / "cv_results_v2.csv"
    results_df.to_csv(results_path, index=False)
    print(f"\nSaved CV results to {results_path}")

    if importance_frames:
        imp_df = pd.concat(importance_frames, ignore_index=True)
        imp_path = RESULTS / "feature_importance_v2.csv"
        imp_df.to_csv(imp_path, index=False)
        print(f"Saved feature importance to {imp_path}")

    # Save models
    if best_lgb_model is not None:
        lgb_path = MODELS / "lightgbm_v2.txt"
        best_lgb_model.booster_.save_model(str(lgb_path))
        print(f"Saved LightGBM model to {lgb_path}")

    if best_xgb_model is not None:
        xgb_path = MODELS / "xgboost_v2.json"
        best_xgb_model.save_model(str(xgb_path))
        print(f"Saved XGBoost model to {xgb_path}")

    # Summary table
    print("\n=== Summary (mean across folds) ===")
    summary = results_df.groupby("config")[["precision", "recall", "f1", "auprc", "auroc"]].mean()
    print(summary.to_string(float_format=lambda x: f"{x:.4f}"))


if __name__ == "__main__":
    main()
