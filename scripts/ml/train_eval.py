"""Train and evaluate ML classifiers on the kam variant calling training dataset.

Runs 5-fold GroupKFold cross-validation with LightGBM and XGBoost,
computes a PASS-filter baseline, and writes results to
docs/project/experiments/02-ml-single-strand/results/.
"""

import sys
from pathlib import Path

import lightgbm as lgb
import numpy as np
import pandas as pd
import xgboost as xgb
from sklearn.metrics import (
    average_precision_score,
    f1_score,
    precision_score,
    recall_score,
    roc_auc_score,
)
from sklearn.model_selection import GroupKFold
from sklearn.preprocessing import LabelEncoder

REPO_ROOT = Path(__file__).resolve().parents[2]
DATA_PATH = REPO_ROOT / "bigdata/experiments/02-ml-single-strand/training_data.csv"
# Small summary outputs (CSV, figures) are committed to docs/.
RESULTS_DIR = REPO_ROOT / "docs/project/experiments/02-ml-single-strand/results"
# Trained model files are large and go to bigdata/.
MODELS_DIR = REPO_ROOT / "bigdata/experiments/02-ml-single-strand/models"

RESULTS_DIR.mkdir(parents=True, exist_ok=True)
MODELS_DIR.mkdir(parents=True, exist_ok=True)

FEATURE_COLS = [
    "vaf", "nref", "nalt", "ndupalt", "nsimalt", "sbp", "conf",
    "duplex_frac", "has_duplex", "ci_width", "ref_len", "alt_len",
    "alt_depth", "variant_class", "mode", "dataset_type",
]
CAT_COLS = ["variant_class", "mode", "dataset_type"]
LABEL_COL = "label"
GROUP_COL = "sample_id"
FILTER_COL = "filter"


def encode_categoricals(df: pd.DataFrame) -> pd.DataFrame:
    """Encode categorical columns with LabelEncoder. Returns a copy."""
    df = df.copy()
    for col in CAT_COLS:
        le = LabelEncoder()
        df[col] = le.fit_transform(df[col].astype(str))
    return df


def baseline_preds(filter_col: pd.Series) -> np.ndarray:
    """Convert filter column to binary predictions (PASS = 1)."""
    return (filter_col == "PASS").astype(int).values


def evaluate_proba(y_true: np.ndarray, y_score: np.ndarray, threshold: float = 0.5) -> dict:
    """Compute precision, recall, F1, AUPRC, AUROC from probability scores."""
    y_pred = (y_score >= threshold).astype(int)
    # Handle edge cases where only one class is present in y_true
    if len(np.unique(y_true)) < 2:
        auprc = float("nan")
        auroc = float("nan")
    else:
        auprc = average_precision_score(y_true, y_score)
        auroc = roc_auc_score(y_true, y_score)
    return {
        "precision": precision_score(y_true, y_pred, zero_division=0),
        "recall": recall_score(y_true, y_pred, zero_division=0),
        "f1": f1_score(y_true, y_pred, zero_division=0),
        "auprc": auprc,
        "auroc": auroc,
    }


def evaluate_binary(y_true: np.ndarray, y_pred: np.ndarray) -> dict:
    """Compute metrics from binary predictions (baseline)."""
    if len(np.unique(y_true)) < 2:
        auprc = float("nan")
        auroc = float("nan")
    else:
        auprc = average_precision_score(y_true, y_pred)
        auroc = roc_auc_score(y_true, y_pred)
    return {
        "precision": precision_score(y_true, y_pred, zero_division=0),
        "recall": recall_score(y_true, y_pred, zero_division=0),
        "f1": f1_score(y_true, y_pred, zero_division=0),
        "auprc": auprc,
        "auroc": auroc,
    }


def main():
    print(f"Loading {DATA_PATH}")
    df = pd.read_csv(DATA_PATH)

    # Drop rows missing any feature or label
    required = FEATURE_COLS + [LABEL_COL, GROUP_COL, FILTER_COL]
    df = df.dropna(subset=required)
    print(f"Rows after dropping missing: {len(df):,}")

    df_enc = encode_categoricals(df)

    X = df_enc[FEATURE_COLS].values
    y = df_enc[LABEL_COL].values
    groups = df_enc[GROUP_COL].values
    filters = df[FILTER_COL].values

    feature_names = FEATURE_COLS

    n_pos = int(y.sum())
    n_neg = int((y == 0).sum())
    scale_pos_weight = n_neg / max(n_pos, 1)
    print(f"Positives: {n_pos}, Negatives: {n_neg}, scale_pos_weight: {scale_pos_weight:.1f}")

    gkf = GroupKFold(n_splits=5)

    records = []
    lgb_models = []
    lgb_auprc_scores = []

    # Store feature importance across folds
    lgb_importance = np.zeros(len(feature_names))

    for fold, (train_idx, val_idx) in enumerate(gkf.split(X, y, groups)):
        X_train, X_val = X[train_idx], X[val_idx]
        y_train, y_val = y[train_idx], y[val_idx]
        filters_val = filters[val_idx]

        fold_n_pos = int(y_train.sum())
        fold_n_neg = int((y_train == 0).sum())
        fold_spw = fold_n_neg / max(fold_n_pos, 1)

        # LightGBM
        lgb_train = lgb.Dataset(X_train, label=y_train, feature_name=feature_names)
        lgb_val = lgb.Dataset(X_val, label=y_val, feature_name=feature_names, reference=lgb_train)

        params_lgb = {
            "objective": "binary",
            "metric": "average_precision",
            "is_unbalance": True,
            "verbose": -1,
            "n_jobs": -1,
            "random_state": 42,
        }
        lgb_model = lgb.train(
            params_lgb,
            lgb_train,
            num_boost_round=200,
            valid_sets=[lgb_val],
            callbacks=[lgb.early_stopping(20, verbose=False), lgb.log_evaluation(-1)],
        )
        lgb_proba = lgb_model.predict(X_val)
        lgb_metrics = evaluate_proba(y_val, lgb_proba)
        lgb_models.append(lgb_model)
        lgb_auprc_scores.append(lgb_metrics["auprc"])

        # Accumulate feature importance (gain)
        imp = lgb_model.feature_importance(importance_type="gain")
        lgb_importance += imp

        # XGBoost
        xgb_train = xgb.DMatrix(X_train, label=y_train, feature_names=feature_names)
        xgb_val = xgb.DMatrix(X_val, label=y_val, feature_names=feature_names)

        params_xgb = {
            "objective": "binary:logistic",
            "eval_metric": "aucpr",
            "scale_pos_weight": fold_spw,
            "verbosity": 0,
            "seed": 42,
            "n_jobs": -1,
        }
        xgb_model = xgb.train(
            params_xgb,
            xgb_train,
            num_boost_round=200,
            evals=[(xgb_val, "val")],
            early_stopping_rounds=20,
            verbose_eval=False,
        )
        xgb_proba = xgb_model.predict(xgb_val)
        xgb_metrics = evaluate_proba(y_val, xgb_proba)

        # Baseline (PASS filter)
        baseline_pred = baseline_preds(pd.Series(filters_val))
        baseline_metrics = evaluate_binary(y_val, baseline_pred)

        for model_name, metrics in [
            ("lightgbm", lgb_metrics),
            ("xgboost", xgb_metrics),
            ("baseline_pass", baseline_metrics),
        ]:
            records.append({
                "model": model_name,
                "fold": fold,
                **metrics,
            })

        print(
            f"Fold {fold} | "
            f"LGB AUPRC={lgb_metrics['auprc']:.4f} F1={lgb_metrics['f1']:.4f} | "
            f"XGB AUPRC={xgb_metrics['auprc']:.4f} F1={xgb_metrics['f1']:.4f} | "
            f"Base AUPRC={baseline_metrics['auprc']:.4f} F1={baseline_metrics['f1']:.4f}"
        )

    # Save CV results
    results_df = pd.DataFrame(records)
    results_path = RESULTS_DIR / "cv_results.csv"
    results_df.to_csv(results_path, index=False)
    print(f"\nCV results written to {results_path}")

    # Save best LightGBM model
    best_fold = int(np.argmax(lgb_auprc_scores))
    best_model = lgb_models[best_fold]
    model_path = MODELS_DIR / "lightgbm_model.txt"
    best_model.save_model(str(model_path))
    print(f"Best LightGBM model (fold {best_fold}, AUPRC={lgb_auprc_scores[best_fold]:.4f}) saved to {model_path}")

    # Feature importance
    avg_importance = lgb_importance / 5.0
    fi_df = pd.DataFrame({
        "feature": feature_names,
        "importance_gain": avg_importance,
    }).sort_values("importance_gain", ascending=False)
    fi_path = RESULTS_DIR / "feature_importance.csv"
    fi_df.to_csv(fi_path, index=False)
    print(f"Feature importance written to {fi_path}")

    # Summary table
    summary = (
        results_df.groupby("model")[["precision", "recall", "f1", "auprc", "auroc"]]
        .mean()
        .round(4)
        .reset_index()
    )
    print("\n=== Cross-validation summary (mean across 5 folds) ===")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
