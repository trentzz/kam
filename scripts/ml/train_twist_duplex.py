#!/usr/bin/env python3
"""Train LightGBM and XGBoost classifiers on the Twist duplex dataset.

Optimises for AUPRC using Optuna HPO (200 trials each).
Outputs calibrated ensemble + ONNX export.

Inputs:
  docs/benchmarking/ml-twist-duplex/train_features.csv.gz
  docs/benchmarking/ml-twist-duplex/test_features.csv.gz

Outputs:
  docs/benchmarking/ml-twist-duplex/models/lightgbm_twist.txt
  docs/benchmarking/ml-twist-duplex/models/xgboost_twist.json
  docs/benchmarking/ml-twist-duplex/models/ensemble_twist.pkl
  docs/benchmarking/ml-twist-duplex/results/cv_results.csv
  docs/benchmarking/ml-twist-duplex/results/test_results.json
  docs/benchmarking/ml-twist-duplex/results/feature_importance_lgb.csv
  docs/benchmarking/ml-twist-duplex/results/feature_importance_xgb.csv
  docs/benchmarking/ml-twist-duplex/results/figures/  (plots)
  docs/benchmarking/ml-twist-duplex/models/xgboost_twist.onnx  (if skl2onnx available)

Usage:
    python3 scripts/ml/train_twist_duplex.py [--hpo-trials N] [--no-onnx]

Run from the repository root.
"""

from __future__ import annotations

import argparse
import json
import pickle
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.calibration import CalibratedClassifierCV
from sklearn.exceptions import UndefinedMetricWarning
from sklearn.metrics import (
    average_precision_score,
    f1_score,
    precision_recall_curve,
    roc_auc_score,
    roc_curve,
)
from sklearn.model_selection import GroupKFold, cross_val_score
from sklearn.preprocessing import LabelEncoder

import lightgbm as lgb
import optuna
import xgboost as xgb

warnings.filterwarnings("ignore", category=UndefinedMetricWarning)
optuna.logging.set_verbosity(optuna.logging.WARNING)

REPO = Path(__file__).resolve().parent.parent.parent
ML_DIR = REPO / "docs" / "benchmarking" / "ml-twist-duplex"
TRAIN_PATH = ML_DIR / "train_features.csv.gz"
TEST_PATH = ML_DIR / "test_features.csv.gz"
RESULTS_DIR = ML_DIR / "results"
FIGURES_DIR = RESULTS_DIR / "figures"
MODELS_DIR = ML_DIR / "models"

RESULTS_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)
MODELS_DIR.mkdir(parents=True, exist_ok=True)

# ─── Feature definitions ──────────────────────────────────────────────────────

# Columns to drop from the feature matrix.
DROP_COLS = {
    "target_id", "ref_seq", "alt_seq", "filter",
    "sample_id", "split", "mode", "label",
    "variant_class_true",  # ground truth type, not a signal
}

# Categorical features for native LightGBM/XGBoost support.
CAT_FEATURES = ["variant_class", "sbp_cat", "depth_bucket"]

VAF_BINS = [
    ("very_low", 0.0, 0.001),
    ("low", 0.001, 0.005),
    ("medium", 0.005, 0.02),
    ("high", 0.02, 1.0),
]

VARIANT_TYPES = [
    "snv", "indel_short", "indel_medium",
    "sv_del", "sv_dup", "sv_inv", "sv_large_del", "ins", "invdel",
]


# ─── Data loading ─────────────────────────────────────────────────────────────

def load_split(path: Path) -> pd.DataFrame:
    """Load a feature CSV and drop rows with ambiguous labels.

    Rows where filter = 'NotTargeted' are removed from tumour-informed data
    because they represent positions outside the target panel and cannot be
    reliably labelled.

    Args:
        path: Path to the CSV.gz file.

    Returns:
        Cleaned DataFrame.
    """
    df = pd.read_csv(path, low_memory=False)
    print(f"Loaded {len(df):,} rows from {path.name}", flush=True)

    # Drop ambiguous filter rows
    if "filter" in df.columns:
        before = len(df)
        df = df[df["filter"] != "NotTargeted"]
        removed = before - len(df)
        if removed:
            print(f"  Dropped {removed:,} 'NotTargeted' rows", flush=True)

    n_pos = (df["label"] == 1).sum()
    n_neg = (df["label"] == 0).sum()
    print(
        f"  Positives: {n_pos:,}  Negatives: {n_neg:,}  "
        f"Positive rate: {n_pos / max(len(df), 1) * 100:.3f}%",
        flush=True,
    )
    return df


def prepare_features(df: pd.DataFrame) -> tuple[pd.DataFrame, list[str]]:
    """Select feature columns and encode categoricals.

    Args:
        df: Raw DataFrame from load_split().

    Returns:
        Tuple (X, feature_cols) where X is the feature matrix and
        feature_cols is the ordered list of column names.
    """
    feature_cols = [c for c in df.columns if c not in DROP_COLS]

    # Encode string categoricals with LabelEncoder
    for col in CAT_FEATURES:
        if col in df.columns and df[col].dtype == object:
            le = LabelEncoder()
            df[col] = le.fit_transform(df[col].astype(str).fillna("unknown"))

    # Fill remaining NaN
    X = df[feature_cols].fillna(0)
    return X, feature_cols


# ─── Metrics ─────────────────────────────────────────────────────────────────

def compute_metrics(
    y_true: np.ndarray,
    y_prob: np.ndarray,
    model_name: str,
    extra: dict | None = None,
) -> dict:
    """Compute a full set of evaluation metrics.

    Args:
        y_true: Ground-truth binary labels.
        y_prob: Predicted probabilities for the positive class.
        model_name: Name for the results dict.
        extra: Optional extra keys to merge.

    Returns:
        Dict of metric values.
    """
    try:
        auroc = float(roc_auc_score(y_true, y_prob))
        auprc = float(average_precision_score(y_true, y_prob))
    except ValueError:
        auroc = auprc = float("nan")

    # Sensitivity at fixed FPR thresholds
    fpr_arr, tpr_arr, thresh_arr = roc_curve(y_true, y_prob)
    sens_at_1pct = float(np.interp(0.01, fpr_arr, tpr_arr))
    sens_at_5pct = float(np.interp(0.05, fpr_arr, tpr_arr))

    # Best F1 threshold
    prec_arr, rec_arr, f1_thresh = precision_recall_curve(y_true, y_prob)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        f1_scores = 2 * prec_arr * rec_arr / (prec_arr + rec_arr + 1e-9)
    best_f1_idx = int(np.argmax(f1_scores))
    best_f1 = float(f1_scores[best_f1_idx])
    best_threshold = float(f1_thresh[min(best_f1_idx, len(f1_thresh) - 1)])

    result = {
        "model": model_name,
        "auroc": auroc,
        "auprc": auprc,
        "sens_at_1pct_fpr": sens_at_1pct,
        "sens_at_5pct_fpr": sens_at_5pct,
        "best_f1": best_f1,
        "best_threshold": best_threshold,
        "n_pos": int(y_true.sum()),
        "n_neg": int((y_true == 0).sum()),
    }
    if extra:
        result.update(extra)
    return result


def per_vaf_bin_metrics(
    df_test: pd.DataFrame,
    y_prob: np.ndarray,
    model_name: str,
) -> list[dict]:
    """Compute AUPRC per VAF bin.

    Args:
        df_test: Test DataFrame with 'vaf_target' and 'label' columns.
        y_prob: Predicted probabilities.
        model_name: Model identifier.

    Returns:
        List of metric dicts, one per VAF bin.
    """
    results = []
    vaf = df_test["vaf_target"].fillna(0).values
    y_true = df_test["label"].values

    for bin_name, vaf_low, vaf_high in VAF_BINS:
        mask = (vaf >= vaf_low) & (vaf < vaf_high)
        if mask.sum() < 10:
            continue
        try:
            auprc = float(average_precision_score(y_true[mask], y_prob[mask]))
        except ValueError:
            auprc = float("nan")
        results.append({
            "model": model_name,
            "vaf_bin": bin_name,
            "vaf_range": f"{vaf_low:.4f}–{vaf_high:.4f}",
            "n_samples": int(mask.sum()),
            "n_pos": int(y_true[mask].sum()),
            "auprc": auprc,
        })
    return results


def per_variant_type_metrics(
    df_test: pd.DataFrame,
    y_prob: np.ndarray,
    model_name: str,
) -> list[dict]:
    """Compute AUPRC per variant type.

    Args:
        df_test: Test DataFrame with 'variant_class_true' and 'label' columns.
        y_prob: Predicted probabilities.
        model_name: Model identifier.

    Returns:
        List of metric dicts, one per variant type.
    """
    results = []
    y_true = df_test["label"].values
    vtypes = df_test["variant_class_true"].fillna("unknown").values

    for vtype in sorted(set(vtypes)):
        mask = vtypes == vtype
        if mask.sum() < 10:
            continue
        try:
            auprc = float(average_precision_score(y_true[mask], y_prob[mask]))
        except ValueError:
            auprc = float("nan")
        results.append({
            "model": model_name,
            "variant_type": vtype,
            "n_samples": int(mask.sum()),
            "n_pos": int(y_true[mask].sum()),
            "auprc": auprc,
        })
    return results


# ─── Optuna HPO: LightGBM ─────────────────────────────────────────────────────

def make_lgb_objective(X_train: np.ndarray, y_train: np.ndarray, groups: np.ndarray):
    """Return an Optuna objective function for LightGBM HPO.

    Uses 3-fold GroupKFold CV for speed during the HPO search.

    Args:
        X_train: Feature matrix.
        y_train: Labels.
        groups: Group identifiers (sample_id) for CV splitting.

    Returns:
        Callable suitable for optuna.Study.optimize().
    """
    def objective(trial: optuna.Trial) -> float:
        """Optimise AUPRC for LightGBM."""
        params = {
            "n_estimators": trial.suggest_int("n_estimators", 1000, 8000),
            "num_leaves": trial.suggest_int("num_leaves", 63, 511),
            "max_depth": trial.suggest_int("max_depth", 6, 15),
            "learning_rate": trial.suggest_float("learning_rate", 0.001, 0.05, log=True),
            "feature_fraction": trial.suggest_float("feature_fraction", 0.5, 1.0),
            "bagging_fraction": trial.suggest_float("bagging_fraction", 0.5, 1.0),
            "bagging_freq": trial.suggest_int("bagging_freq", 1, 10),
            "lambda_l1": trial.suggest_float("lambda_l1", 0.0, 10.0),
            "lambda_l2": trial.suggest_float("lambda_l2", 0.0, 10.0),
            "min_child_samples": trial.suggest_int("min_child_samples", 10, 100),
            "min_child_weight": trial.suggest_float("min_child_weight", 0.001, 10.0, log=True),
            "is_unbalance": True,
            "objective": "binary",
            "metric": "average_precision",
            "verbosity": -1,
            "n_jobs": -1,
            "random_state": 42,
        }
        model = lgb.LGBMClassifier(**params)
        scores = cross_val_score(
            model, X_train, y_train,
            groups=groups,
            cv=GroupKFold(n_splits=3),
            scoring="average_precision",
        )
        return float(scores.mean())

    return objective


# ─── Optuna HPO: XGBoost ──────────────────────────────────────────────────────

def make_xgb_objective(
    X_train: np.ndarray,
    y_train: np.ndarray,
    groups: np.ndarray,
    n_pos: int,
    n_neg: int,
):
    """Return an Optuna objective function for XGBoost HPO.

    Args:
        X_train: Feature matrix.
        y_train: Labels.
        groups: Group identifiers.
        n_pos: Number of positive training examples.
        n_neg: Number of negative training examples.

    Returns:
        Callable for optuna.Study.optimize().
    """
    scale_pos_weight = n_neg / max(n_pos, 1)

    def objective(trial: optuna.Trial) -> float:
        """Optimise AUPRC for XGBoost."""
        params = {
            "n_estimators": trial.suggest_int("n_estimators", 1000, 8000),
            "max_depth": trial.suggest_int("max_depth", 6, 12),
            "learning_rate": trial.suggest_float("learning_rate", 0.001, 0.05, log=True),
            "subsample": trial.suggest_float("subsample", 0.5, 1.0),
            "colsample_bytree": trial.suggest_float("colsample_bytree", 0.5, 1.0),
            "colsample_bylevel": trial.suggest_float("colsample_bylevel", 0.5, 1.0),
            "min_child_weight": trial.suggest_int("min_child_weight", 1, 50),
            "gamma": trial.suggest_float("gamma", 0.0, 5.0),
            "reg_alpha": trial.suggest_float("reg_alpha", 0.0, 10.0),
            "reg_lambda": trial.suggest_float("reg_lambda", 0.0, 10.0),
            "scale_pos_weight": scale_pos_weight,
            "tree_method": "hist",
            "eval_metric": "aucpr",
            "use_label_encoder": False,
            "verbosity": 0,
            "n_jobs": -1,
            "random_state": 42,
        }
        model = xgb.XGBClassifier(**params)
        scores = cross_val_score(
            model, X_train, y_train,
            groups=groups,
            cv=GroupKFold(n_splits=3),
            scoring="average_precision",
        )
        return float(scores.mean())

    return objective


# ─── 5-fold CV ────────────────────────────────────────────────────────────────

def run_cv_final(
    model,
    X: np.ndarray,
    y: np.ndarray,
    groups: np.ndarray,
    model_name: str,
    n_splits: int = 5,
) -> list[dict]:
    """Run 5-fold GroupKFold CV and return per-fold metrics.

    Args:
        model: Scikit-learn compatible model.
        X: Feature matrix.
        y: Labels.
        groups: Group identifiers.
        model_name: Name for the results dict.
        n_splits: Number of CV folds.

    Returns:
        List of metric dicts, one per fold.
    """
    gkf = GroupKFold(n_splits=n_splits)
    fold_results = []
    for fold, (train_idx, val_idx) in enumerate(gkf.split(X, y, groups)):
        model.fit(X[train_idx], y[train_idx])
        y_prob = model.predict_proba(X[val_idx])[:, 1]
        m = compute_metrics(y[val_idx], y_prob, model_name)
        m["fold"] = fold
        fold_results.append(m)
        print(
            f"  fold {fold}: AUROC={m['auroc']:.4f} AUPRC={m['auprc']:.4f} "
            f"F1={m['best_f1']:.4f}",
            flush=True,
        )
    return fold_results


# ─── Plotting ─────────────────────────────────────────────────────────────────

def save_plots(
    y_test: np.ndarray,
    probs_dict: dict[str, np.ndarray],
    df_test: pd.DataFrame,
) -> None:
    """Save PR curves, ROC curves, and VAF-bin breakdown plots.

    Skips gracefully if matplotlib is not available.

    Args:
        y_test: Ground-truth test labels.
        probs_dict: Dict of model_name -> probability array.
        df_test: Test DataFrame with vaf_target column.
    """
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("[WARN] matplotlib not available — skipping plots.", file=sys.stderr)
        return

    # PR curves
    fig, ax = plt.subplots(figsize=(7, 6))
    for name, probs in probs_dict.items():
        prec, rec, _ = precision_recall_curve(y_test, probs)
        ap = average_precision_score(y_test, probs)
        ax.plot(rec, prec, label=f"{name} (AUPRC={ap:.3f})")
    ax.set_xlabel("Recall")
    ax.set_ylabel("Precision")
    ax.set_title("Precision-Recall Curve — Twist Duplex Test Set")
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "pr_curves.png", dpi=150)
    plt.close(fig)

    # ROC curves
    fig, ax = plt.subplots(figsize=(7, 6))
    for name, probs in probs_dict.items():
        fpr, tpr, _ = roc_curve(y_test, probs)
        auroc = roc_auc_score(y_test, probs)
        ax.plot(fpr, tpr, label=f"{name} (AUROC={auroc:.3f})")
    ax.plot([0, 1], [0, 1], "k--", alpha=0.4)
    ax.set_xlabel("False Positive Rate")
    ax.set_ylabel("True Positive Rate")
    ax.set_title("ROC Curve — Twist Duplex Test Set")
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "roc_curves.png", dpi=150)
    plt.close(fig)

    # AUPRC by VAF bin (ensemble)
    if "ensemble" in probs_dict:
        bins = []
        vaf = df_test["vaf_target"].fillna(0).values
        probs = probs_dict["ensemble"]
        y = y_test
        for bin_name, vaf_low, vaf_high in VAF_BINS:
            mask = (vaf >= vaf_low) & (vaf < vaf_high)
            if mask.sum() < 10:
                continue
            try:
                ap = average_precision_score(y[mask], probs[mask])
            except ValueError:
                ap = 0.0
            bins.append((bin_name, ap, mask.sum()))

        if bins:
            fig, ax = plt.subplots(figsize=(8, 5))
            names = [b[0] for b in bins]
            values = [b[1] for b in bins]
            counts = [b[2] for b in bins]
            bars = ax.bar(names, values, color="steelblue", alpha=0.8)
            for bar, cnt in zip(bars, counts):
                ax.text(
                    bar.get_x() + bar.get_width() / 2,
                    bar.get_height() + 0.01,
                    f"n={cnt}",
                    ha="center", va="bottom", fontsize=9,
                )
            ax.set_xlabel("VAF Bin")
            ax.set_ylabel("AUPRC")
            ax.set_title("AUPRC by VAF Bin (Ensemble) — Twist Duplex Test Set")
            ax.set_ylim(0, 1.05)
            ax.grid(True, axis="y", alpha=0.3)
            fig.tight_layout()
            fig.savefig(FIGURES_DIR / "auprc_by_vaf.png", dpi=150)
            plt.close(fig)

    print(f"Figures saved to {FIGURES_DIR}", flush=True)


# ─── ONNX export ──────────────────────────────────────────────────────────────

def export_onnx(model, feature_cols: list[str]) -> None:
    """Export an XGBoost model to ONNX format.

    Requires skl2onnx. Skips gracefully if not installed.

    Args:
        model: Trained XGBoost classifier.
        feature_cols: Ordered list of feature column names.
    """
    try:
        from skl2onnx import convert_sklearn
        from skl2onnx.common.data_types import FloatTensorType
    except ImportError:
        print("[WARN] skl2onnx not available — skipping ONNX export.", file=sys.stderr)
        return

    initial_type = [("float_input", FloatTensorType([None, len(feature_cols)]))]
    try:
        onnx_model = convert_sklearn(model, initial_types=initial_type)
        out_path = MODELS_DIR / "xgboost_twist.onnx"
        with open(out_path, "wb") as fh:
            fh.write(onnx_model.SerializeToString())
        print(f"ONNX model written: {out_path}", flush=True)
    except Exception as exc:
        print(f"[WARN] ONNX export failed: {exc}", file=sys.stderr)


# ─── CLI ──────────────────────────────────────────────────────────────────────

def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Train LightGBM and XGBoost classifiers on the Twist duplex dataset.",
    )
    parser.add_argument(
        "--hpo-trials",
        type=int,
        default=200,
        help="Number of Optuna HPO trials per model (default: 200).",
    )
    parser.add_argument(
        "--no-onnx",
        action="store_true",
        help="Skip ONNX export.",
    )
    parser.add_argument(
        "--n-jobs",
        type=int,
        default=-1,
        help="Parallel jobs for model training (default: -1 = all cores).",
    )
    return parser.parse_args()


# ─── Main ─────────────────────────────────────────────────────────────────────

def main() -> None:
    """Entry point."""
    args = parse_args()

    for path in (TRAIN_PATH, TEST_PATH):
        if not path.exists():
            print(
                f"[ERROR] {path} not found.\n"
                "Run scripts/ml/build_twist_duplex_features.py first.",
                file=sys.stderr,
            )
            sys.exit(1)

    # ── Load data ─────────────────────────────────────────────────────────────
    print("\n=== Loading data ===", flush=True)
    df_train = load_split(TRAIN_PATH)
    df_test = load_split(TEST_PATH)

    X_train_df, feature_cols = prepare_features(df_train)
    X_test_df, _ = prepare_features(df_test)

    X_train = X_train_df.values
    y_train = df_train["label"].values
    groups = df_train["sample_id"].values if "sample_id" in df_train.columns else np.arange(len(y_train))

    X_test = X_test_df.values
    y_test = df_test["label"].values

    n_pos = int(y_train.sum())
    n_neg = int((y_train == 0).sum())
    print(f"Train: {len(X_train):,} rows, {n_pos:,} pos, {n_neg:,} neg", flush=True)
    print(f"Test:  {len(X_test):,} rows, {int(y_test.sum()):,} pos", flush=True)
    print(f"Features: {len(feature_cols)}", flush=True)

    all_cv_results: list[dict] = []
    test_results: dict[str, dict] = {}

    # ── LightGBM HPO ──────────────────────────────────────────────────────────
    print(f"\n=== LightGBM HPO ({args.hpo_trials} trials) ===", flush=True)
    lgb_study = optuna.create_study(direction="maximize", sampler=optuna.samplers.TPESampler(seed=42))
    lgb_study.optimize(
        make_lgb_objective(X_train, y_train, groups),
        n_trials=args.hpo_trials,
        show_progress_bar=True,
    )
    lgb_best_params = lgb_study.best_params
    lgb_best_params.update({
        "is_unbalance": True,
        "objective": "binary",
        "metric": "average_precision",
        "verbosity": -1,
        "n_jobs": args.n_jobs,
        "random_state": 42,
    })
    print(f"Best LightGBM AUPRC (CV): {lgb_study.best_value:.4f}", flush=True)

    # ── XGBoost HPO ───────────────────────────────────────────────────────────
    print(f"\n=== XGBoost HPO ({args.hpo_trials} trials) ===", flush=True)
    xgb_study = optuna.create_study(direction="maximize", sampler=optuna.samplers.TPESampler(seed=42))
    xgb_study.optimize(
        make_xgb_objective(X_train, y_train, groups, n_pos, n_neg),
        n_trials=args.hpo_trials,
        show_progress_bar=True,
    )
    xgb_best_params = xgb_study.best_params
    xgb_best_params.update({
        "scale_pos_weight": n_neg / max(n_pos, 1),
        "tree_method": "hist",
        "eval_metric": "aucpr",
        "use_label_encoder": False,
        "verbosity": 0,
        "n_jobs": args.n_jobs,
        "random_state": 42,
    })
    print(f"Best XGBoost AUPRC (CV): {xgb_study.best_value:.4f}", flush=True)

    # ── 5-fold CV on full train set ───────────────────────────────────────────
    print("\n=== 5-fold GroupKFold CV ===", flush=True)

    print("\nLightGBM:", flush=True)
    lgb_model_cv = lgb.LGBMClassifier(**lgb_best_params)
    lgb_cv_results = run_cv_final(lgb_model_cv, X_train, y_train, groups, "lightgbm")
    all_cv_results.extend(lgb_cv_results)

    print("\nXGBoost:", flush=True)
    xgb_model_cv = xgb.XGBClassifier(**xgb_best_params)
    xgb_cv_results = run_cv_final(xgb_model_cv, X_train, y_train, groups, "xgboost")
    all_cv_results.extend(xgb_cv_results)

    cv_df = pd.DataFrame(all_cv_results)
    cv_df.to_csv(RESULTS_DIR / "cv_results.csv", index=False)
    print(f"\nCV results written: {RESULTS_DIR / 'cv_results.csv'}", flush=True)

    # ── Final training on full train set ──────────────────────────────────────
    print("\n=== Final training on full training set ===", flush=True)

    lgb_model = lgb.LGBMClassifier(**lgb_best_params)
    lgb_model.fit(X_train, y_train)
    lgb_model.booster_.save_model(str(MODELS_DIR / "lightgbm_twist.txt"))
    print(f"LightGBM model saved: {MODELS_DIR / 'lightgbm_twist.txt'}", flush=True)

    xgb_model = xgb.XGBClassifier(**xgb_best_params)
    xgb_model.fit(X_train, y_train)
    xgb_model.save_model(str(MODELS_DIR / "xgboost_twist.json"))
    print(f"XGBoost model saved: {MODELS_DIR / 'xgboost_twist.json'}", flush=True)

    # ── Calibration ───────────────────────────────────────────────────────────
    print("\n=== Calibrating models ===", flush=True)
    lgb_cal = CalibratedClassifierCV(lgb_model, method="isotonic", cv="prefit")
    lgb_cal.fit(X_train, y_train)

    xgb_cal = CalibratedClassifierCV(xgb_model, method="isotonic", cv="prefit")
    xgb_cal.fit(X_train, y_train)

    # ── Test set evaluation ───────────────────────────────────────────────────
    print("\n=== Test set evaluation ===", flush=True)
    prob_lgb = lgb_cal.predict_proba(X_test)[:, 1]
    prob_xgb = xgb_cal.predict_proba(X_test)[:, 1]
    prob_ensemble = 0.5 * prob_lgb + 0.5 * prob_xgb

    probs_dict = {
        "lightgbm": prob_lgb,
        "xgboost": prob_xgb,
        "ensemble": prob_ensemble,
    }

    all_test_results: list[dict] = []
    for model_name, probs in probs_dict.items():
        m = compute_metrics(y_test, probs, model_name)
        test_results[model_name] = m
        print(
            f"  {model_name:<12} AUROC={m['auroc']:.4f} AUPRC={m['auprc']:.4f} "
            f"F1={m['best_f1']:.4f}",
            flush=True,
        )
        all_test_results.append(m)
        all_test_results.extend(per_vaf_bin_metrics(df_test, probs, model_name))
        all_test_results.extend(per_variant_type_metrics(df_test, probs, model_name))

    with open(RESULTS_DIR / "test_results.json", "w") as fh:
        json.dump(all_test_results, fh, indent=2)
    print(f"\nTest results written: {RESULTS_DIR / 'test_results.json'}", flush=True)

    # ── Feature importance ────────────────────────────────────────────────────
    lgb_imp = pd.DataFrame({
        "feature": feature_cols,
        "importance": lgb_model.feature_importances_,
    }).sort_values("importance", ascending=False)
    lgb_imp.to_csv(RESULTS_DIR / "feature_importance_lgb.csv", index=False)

    xgb_imp = pd.DataFrame({
        "feature": feature_cols,
        "importance": xgb_model.feature_importances_,
    }).sort_values("importance", ascending=False)
    xgb_imp.to_csv(RESULTS_DIR / "feature_importance_xgb.csv", index=False)

    print(f"\nTop 10 LightGBM features:")
    for _, row in lgb_imp.head(10).iterrows():
        print(f"  {row['feature']:<40} {row['importance']:.0f}")

    # ── Save ensemble ─────────────────────────────────────────────────────────
    ensemble = {"lgb": lgb_cal, "xgb": xgb_cal, "feature_cols": feature_cols}
    with open(MODELS_DIR / "ensemble_twist.pkl", "wb") as fh:
        pickle.dump(ensemble, fh)
    print(f"\nEnsemble saved: {MODELS_DIR / 'ensemble_twist.pkl'}", flush=True)

    # ── ONNX export ───────────────────────────────────────────────────────────
    if not args.no_onnx:
        print("\n=== ONNX export ===", flush=True)
        export_onnx(xgb_model, feature_cols)

    # ── Save plots ────────────────────────────────────────────────────────────
    print("\n=== Saving plots ===", flush=True)
    save_plots(y_test, probs_dict, df_test)

    print("\n=== Training complete ===", flush=True)


if __name__ == "__main__":
    main()
