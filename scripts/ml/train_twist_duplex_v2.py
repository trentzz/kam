#!/usr/bin/env python3
"""Train LightGBM and XGBoost classifiers for the twist-duplex-v2 model.

Uses HPO via Optuna, 4-fold GroupKFold CV, isotonic calibration, and ONNX
export. Produces model metadata JSON compatible with the kam Rust inference
pipeline.

Inputs:
  input_dir/train_features.csv.gz
  input_dir/test_features.csv.gz

Outputs:
  output_dir/models/twist-duplex-v2.onnx
  output_dir/models/twist-duplex-v2.json
  output_dir/results/cv_results_v2.csv
  output_dir/results/test_results_v2.json
  output_dir/results/feature_importance_lgb_v2.csv
  results_dir/pr_curve_v2.png
  results_dir/ml_prob_histogram_v2.png
  results_dir/feature_importance_v2.png

Usage:
    python3 scripts/ml/train_twist_duplex_v2.py \\
        --input-dir bigdata/experiments/03-ml-twist-duplex-v2/ \\
        --output-dir bigdata/experiments/03-ml-twist-duplex-v2/ \\
        --results-dir docs/project/experiments/03-ml-twist-duplex-v2/results/ \\
        --hpo-trials 50

Run from the repository root.
"""

from __future__ import annotations

import argparse
import json
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
)
from sklearn.model_selection import GroupKFold, cross_val_score

import lightgbm as lgb
import optuna
import xgboost as xgb

warnings.filterwarnings("ignore", category=UndefinedMetricWarning)
optuna.logging.set_verbosity(optuna.logging.WARNING)

# ── Feature list: must match build_real_training_data.py exactly ──────────────

FEATURE_NAMES: list[str] = [
    # Original 33
    "vaf",
    "nref",
    "nalt",
    "ndupalt",
    "nsimalt",
    "sbp",
    "conf",
    "ref_len",
    "alt_len",
    "duplex_frac",
    "has_duplex",
    "ci_width",
    "alt_depth",
    "log_nalt",
    "log_nref",
    "log_alt_depth",
    "log_vaf",
    "vaf_times_conf",
    "vaf_times_nalt",
    "nalt_over_conf",
    "ci_width_rel",
    "snr",
    "conf_sq",
    "nalt_sq",
    "vaf_sq",
    "ref_alt_len_ratio",
    "indel_size",
    "duplex_enrichment",
    "simplex_only_frac",
    "conf_above_99",
    "conf_above_999",
    "sbp_above_05",
    "variant_class_enc",
    # Category B
    "n_simplex_fwd_alt",
    "n_simplex_rev_alt",
    "n_duplex_ref",
    "n_simplex_ref",
    "mean_alt_error_prob",
    "min_variant_specific_duplex",
    "mean_variant_specific_molecules",
    # Category C
    "strand_asymmetry_alt",
    "duplex_vaf",
    "simplex_vaf",
    "duplex_simplex_vaf_delta",
    # Category A
    "subst_type",
    "trinuc_context",
    "is_cpg",
    "gc_content_ref",
    "homopolymer_run",
]

# Metadata columns to exclude from the feature matrix.
META_COLS = {"label", "sample_id", "ng_condition", "vaf_nominal"}

# VAF bins for per-bin AUPRC reporting.
VAF_BINS: list[tuple[str, float, float]] = [
    ("very_low", 0.0, 0.001),
    ("low", 0.001, 0.005),
    ("medium", 0.005, 0.02),
    ("high", 0.02, 1.0),
]

VARIANT_CLASS_MAP: dict[str, int] = {
    "SNV": 0,
    "Insertion": 1,
    "Deletion": 2,
    "MNV": 3,
    "Complex": 4,
}


# ── Data loading ───────────────────────────────────────────────────────────────

def load_split(path: Path) -> pd.DataFrame:
    """Load a compressed feature CSV.

    Args:
        path: Path to the csv.gz file.

    Returns:
        DataFrame with all feature and metadata columns.

    Example:
        >>> df = load_split(Path("train_features.csv.gz"))
        >>> "label" in df.columns
        True
    """
    df = pd.read_csv(path, low_memory=False)
    print(f"Loaded {len(df):,} rows from {path.name}", flush=True)
    n_pos = int((df["label"] == 1).sum())
    n_neg = int((df["label"] == 0).sum())
    print(
        f"  Positives: {n_pos:,}  Negatives: {n_neg:,}  "
        f"Positive rate: {n_pos / max(len(df), 1) * 100:.3f}%",
        flush=True,
    )
    return df


def extract_arrays(
    df: pd.DataFrame,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Extract feature matrix, label vector, and group vector from a DataFrame.

    Fills NaN values with 0. The group vector uses sample_id for GroupKFold.

    Args:
        df: DataFrame from load_split().

    Returns:
        Tuple (X, y, groups) where X has shape (n, 49), y has shape (n,),
        and groups has shape (n,).

    Example:
        >>> X, y, g = extract_arrays(df)
        >>> X.shape[1] == 49
        True
    """
    X = df[FEATURE_NAMES].fillna(0).values
    y = df["label"].values.astype(int)
    groups = df["sample_id"].values if "sample_id" in df.columns else np.arange(len(y))
    return X, y, groups


# ── Metrics ────────────────────────────────────────────────────────────────────

def compute_threshold_metrics(
    y_true: np.ndarray,
    y_prob: np.ndarray,
    thresholds: list[float],
) -> list[dict]:
    """Compute precision, recall, and FP count at fixed thresholds.

    Args:
        y_true: Ground-truth labels.
        y_prob: Predicted probabilities.
        thresholds: List of probability thresholds to evaluate.

    Returns:
        List of dicts with keys: threshold, precision, recall, fp_count, tp_count.

    Example:
        >>> rows = compute_threshold_metrics(y_true, y_prob, [0.5])
        >>> "fp_count" in rows[0]
        True
    """
    results = []
    for thresh in thresholds:
        pred = (y_prob >= thresh).astype(int)
        tp = int(((pred == 1) & (y_true == 1)).sum())
        fp = int(((pred == 1) & (y_true == 0)).sum())
        fn = int(((pred == 0) & (y_true == 1)).sum())
        prec = tp / (tp + fp) if (tp + fp) > 0 else 0.0
        rec = tp / (tp + fn) if (tp + fn) > 0 else 0.0
        results.append({
            "threshold": thresh,
            "precision": prec,
            "recall": rec,
            "tp_count": tp,
            "fp_count": fp,
        })
    return results


def find_best_f1_threshold(y_true: np.ndarray, y_prob: np.ndarray) -> tuple[float, float]:
    """Find the probability threshold that maximises F1 score.

    Args:
        y_true: Ground-truth labels.
        y_prob: Predicted probabilities.

    Returns:
        Tuple (best_threshold, best_f1).

    Example:
        >>> thresh, f1 = find_best_f1_threshold(y_true, y_prob)
        >>> 0.0 <= thresh <= 1.0
        True
    """
    prec_arr, rec_arr, thresh_arr = precision_recall_curve(y_true, y_prob)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        f1_scores = 2 * prec_arr * rec_arr / (prec_arr + rec_arr + 1e-9)
    best_idx = int(np.argmax(f1_scores))
    best_f1 = float(f1_scores[best_idx])
    best_threshold = float(thresh_arr[min(best_idx, len(thresh_arr) - 1)])
    return best_threshold, best_f1


def per_vaf_bin_auprc(
    df_test: pd.DataFrame,
    y_prob: np.ndarray,
    model_name: str,
) -> list[dict]:
    """Compute AUPRC within each nominal VAF bin.

    Args:
        df_test: Test DataFrame with 'vaf_nominal' and 'label' columns.
        y_prob: Predicted probabilities.
        model_name: Label for the results dict.

    Returns:
        List of dicts with keys: model, vaf_bin, vaf_range, n_samples, n_pos, auprc.

    Example:
        >>> rows = per_vaf_bin_auprc(df_test, probs, "lgb")
        >>> all("auprc" in r for r in rows)
        True
    """
    results = []
    vaf = df_test["vaf_nominal"].fillna(0).values
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
            "vaf_range": f"{vaf_low:.4f}-{vaf_high:.4f}",
            "n_samples": int(mask.sum()),
            "n_pos": int(y_true[mask].sum()),
            "auprc": auprc,
        })
    return results


# ── HPO objectives ─────────────────────────────────────────────────────────────

def make_lgb_objective(
    X_train: np.ndarray,
    y_train: np.ndarray,
    groups: np.ndarray,
):
    """Return an Optuna objective for LightGBM AUPRC optimisation.

    Uses 4-fold GroupKFold during HPO for speed.

    Args:
        X_train: Feature matrix.
        y_train: Labels.
        groups: Sample group identifiers for CV splitting.

    Returns:
        Callable for optuna.Study.optimize().

    Example:
        >>> obj = make_lgb_objective(X, y, g)
        >>> callable(obj)
        True
    """
    def objective(trial: optuna.Trial) -> float:
        """Optimise mean AUPRC across 4 GroupKFold folds."""
        params = {
            "n_estimators": trial.suggest_int("n_estimators", 500, 5000),
            "num_leaves": trial.suggest_int("num_leaves", 31, 255),
            "max_depth": trial.suggest_int("max_depth", 4, 12),
            "learning_rate": trial.suggest_float("learning_rate", 0.005, 0.1, log=True),
            "feature_fraction": trial.suggest_float("feature_fraction", 0.6, 1.0),
            "bagging_fraction": trial.suggest_float("bagging_fraction", 0.6, 1.0),
            "bagging_freq": 5,
            "min_child_samples": trial.suggest_int("min_child_samples", 10, 100),
            "reg_alpha": trial.suggest_float("reg_alpha", 1e-8, 10.0, log=True),
            "reg_lambda": trial.suggest_float("reg_lambda", 1e-8, 10.0, log=True),
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
            cv=GroupKFold(n_splits=4),
            scoring="average_precision",
        )
        return float(scores.mean())

    return objective


def make_xgb_objective(
    X_train: np.ndarray,
    y_train: np.ndarray,
    groups: np.ndarray,
    scale_pos_weight: float,
):
    """Return an Optuna objective for XGBoost AUPRC optimisation.

    Args:
        X_train: Feature matrix.
        y_train: Labels.
        groups: Sample group identifiers for CV splitting.
        scale_pos_weight: n_neg / n_pos to address class imbalance.

    Returns:
        Callable for optuna.Study.optimize().

    Example:
        >>> obj = make_xgb_objective(X, y, g, 10.0)
        >>> callable(obj)
        True
    """
    def objective(trial: optuna.Trial) -> float:
        """Optimise mean AUPRC across 4 GroupKFold folds."""
        params = {
            "n_estimators": trial.suggest_int("n_estimators", 500, 5000),
            "max_depth": trial.suggest_int("max_depth", 4, 12),
            "learning_rate": trial.suggest_float("learning_rate", 0.005, 0.1, log=True),
            "subsample": trial.suggest_float("subsample", 0.6, 1.0),
            "colsample_bytree": trial.suggest_float("colsample_bytree", 0.6, 1.0),
            "min_child_weight": trial.suggest_int("min_child_weight", 10, 100),
            "reg_alpha": trial.suggest_float("reg_alpha", 1e-8, 10.0, log=True),
            "reg_lambda": trial.suggest_float("reg_lambda", 1e-8, 10.0, log=True),
            "scale_pos_weight": scale_pos_weight,
            "tree_method": "hist",
            "eval_metric": "aucpr",
            "verbosity": 0,
            "n_jobs": -1,
            "random_state": 42,
        }
        model = xgb.XGBClassifier(**params)
        scores = cross_val_score(
            model, X_train, y_train,
            groups=groups,
            cv=GroupKFold(n_splits=4),
            scoring="average_precision",
        )
        return float(scores.mean())

    return objective


# ── Cross-validation ───────────────────────────────────────────────────────────

def run_cv(
    model,
    X: np.ndarray,
    y: np.ndarray,
    groups: np.ndarray,
    model_name: str,
    n_splits: int = 4,
) -> list[dict]:
    """Run 4-fold GroupKFold CV and return per-fold metrics.

    Args:
        model: Scikit-learn compatible classifier.
        X: Feature matrix.
        y: Labels.
        groups: Group identifiers for CV splitting.
        model_name: Label for results dicts.
        n_splits: Number of folds.

    Returns:
        List of dicts with fold, model, auroc, auprc, best_f1.

    Example:
        >>> results = run_cv(lgb_model, X, y, groups, "lgb")
        >>> len(results) == 4
        True
    """
    gkf = GroupKFold(n_splits=n_splits)
    fold_results = []
    for fold, (train_idx, val_idx) in enumerate(gkf.split(X, y, groups)):
        model.fit(X[train_idx], y[train_idx])
        y_prob = model.predict_proba(X[val_idx])[:, 1]
        try:
            auroc = float(roc_auc_score(y[val_idx], y_prob))
            auprc = float(average_precision_score(y[val_idx], y_prob))
        except ValueError:
            auroc = auprc = float("nan")
        thresh, best_f1 = find_best_f1_threshold(y[val_idx], y_prob)
        row = {
            "fold": fold,
            "model": model_name,
            "auroc": auroc,
            "auprc": auprc,
            "best_f1": best_f1,
            "best_threshold": thresh,
        }
        fold_results.append(row)
        print(
            f"  fold {fold}: AUROC={auroc:.4f}  AUPRC={auprc:.4f}  F1={best_f1:.4f}",
            flush=True,
        )
    return fold_results


# ── ONNX export ────────────────────────────────────────────────────────────────

def export_onnx_lgb(model, n_features: int, out_path: Path) -> bool:
    """Export a LightGBM classifier to ONNX via skl2onnx.

    Falls back to XGBoost export if skl2onnx is not available.

    Args:
        model: Trained LightGBM sklearn classifier or calibrated wrapper.
        n_features: Number of input features.
        out_path: Destination path for the .onnx file.

    Returns:
        True if export succeeded, False otherwise.

    Example:
        >>> ok = export_onnx_lgb(lgb_model, 49, Path("model.onnx"))
        >>> isinstance(ok, bool)
        True
    """
    try:
        import onnxmltools
        from onnxmltools.convert.common.data_types import FloatTensorType as OnnxFloat
    except ImportError:
        print("[WARN] onnxmltools not available; skipping ONNX export.", file=sys.stderr)
        return False

    # Export the raw booster, not the calibrated wrapper.
    # CalibratedClassifierCV wraps the estimator; unwrap it.
    raw = model
    if hasattr(raw, "estimator"):
        raw = raw.estimator
    if hasattr(raw, "base_estimator"):
        raw = raw.base_estimator
    booster = raw.booster_ if hasattr(raw, "booster_") else raw

    try:
        onnx_model = onnxmltools.convert_lightgbm(
            booster,
            initial_types=[("float_input", OnnxFloat([None, n_features]))],
        )
        out_path.parent.mkdir(parents=True, exist_ok=True)
        onnxmltools.utils.save_model(onnx_model, str(out_path))
        print(f"ONNX model written: {out_path}", flush=True)
        return True
    except Exception as exc:
        print(f"[WARN] ONNX export failed: {exc}", file=sys.stderr)
        return False


def export_onnx_xgb(model, n_features: int, out_path: Path) -> bool:
    """Export an XGBoost classifier to ONNX via skl2onnx (fallback).

    Args:
        model: Trained XGBoost sklearn classifier or calibrated wrapper.
        n_features: Number of input features.
        out_path: Destination path for the .onnx file.

    Returns:
        True if export succeeded, False otherwise.

    Example:
        >>> ok = export_onnx_xgb(xgb_model, 49, Path("model.onnx"))
        >>> isinstance(ok, bool)
        True
    """
    try:
        from skl2onnx import convert_sklearn
        from skl2onnx.common.data_types import FloatTensorType
    except ImportError:
        print("[WARN] skl2onnx not available; skipping XGBoost ONNX export.", file=sys.stderr)
        return False

    initial_type = [("float_input", FloatTensorType([None, n_features]))]
    try:
        onnx_model = convert_sklearn(model, initial_types=initial_type)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        with open(out_path, "wb") as fh:
            fh.write(onnx_model.SerializeToString())
        print(f"ONNX model (XGBoost fallback) written: {out_path}", flush=True)
        return True
    except Exception as exc:
        print(f"[WARN] XGBoost ONNX export failed: {exc}", file=sys.stderr)
        return False


# ── Plotting ───────────────────────────────────────────────────────────────────

def save_plots(
    y_test: np.ndarray,
    prob_lgb: np.ndarray,
    results_dir: Path,
    feature_imp_df: pd.DataFrame,
) -> None:
    """Save probability histogram, PR curve, and feature importance bar chart.

    Skips gracefully if matplotlib is unavailable.

    Args:
        y_test: Ground-truth test labels.
        prob_lgb: LightGBM predicted probabilities on the test set.
        results_dir: Directory to write figures.
        feature_imp_df: DataFrame with 'feature' and 'importance' columns.

    Example:
        >>> save_plots(y_test, probs, Path("results/"), imp_df)
    """
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("[WARN] matplotlib not available; skipping plots.", file=sys.stderr)
        return

    results_dir.mkdir(parents=True, exist_ok=True)

    # ── ml_prob histogram ────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(8, 5))
    tp_mask = y_test == 1
    fp_mask = y_test == 0
    bins = np.linspace(0, 1, 41)
    ax.hist(prob_lgb[tp_mask], bins=bins, alpha=0.7, label="TP", color="steelblue")
    ax.hist(prob_lgb[fp_mask], bins=bins, alpha=0.7, label="FP", color="tomato")
    ax.set_xlabel("ml_prob")
    ax.set_ylabel("Count")
    ax.set_title("ML probability distribution — twist-duplex-v2 test set")
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(results_dir / "ml_prob_histogram_v2.png", dpi=150)
    plt.close(fig)

    # ── PR curve ─────────────────────────────────────────────────────────────
    prec, rec, _ = precision_recall_curve(y_test, prob_lgb)
    auprc = average_precision_score(y_test, prob_lgb)
    fig, ax = plt.subplots(figsize=(7, 6))
    ax.plot(rec, prec, label=f"LightGBM (AUPRC={auprc:.3f})", color="steelblue")
    ax.set_xlabel("Recall")
    ax.set_ylabel("Precision")
    ax.set_title("Precision-Recall Curve — twist-duplex-v2 test set")
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(results_dir / "pr_curve_v2.png", dpi=150)
    plt.close(fig)

    # ── Feature importance bar chart (top 20) ─────────────────────────────────
    top20 = feature_imp_df.head(20)
    fig, ax = plt.subplots(figsize=(10, 7))
    ax.barh(top20["feature"][::-1], top20["importance"][::-1], color="steelblue", alpha=0.85)
    ax.set_xlabel("Importance")
    ax.set_title("Top 20 feature importances — LightGBM twist-duplex-v2")
    ax.grid(True, axis="x", alpha=0.3)
    fig.tight_layout()
    fig.savefig(results_dir / "feature_importance_v2.png", dpi=150)
    plt.close(fig)

    print(f"Figures written to {results_dir}", flush=True)


# ── CLI ────────────────────────────────────────────────────────────────────────

def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="HPO + train + calibrate + ONNX export for twist-duplex-v2.",
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        required=True,
        help="Directory containing train_features.csv.gz and test_features.csv.gz.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Directory for models/ and results/ output.",
    )
    parser.add_argument(
        "--results-dir",
        type=Path,
        required=True,
        help="Directory for figures (PR curve, histogram, feature importance).",
    )
    parser.add_argument(
        "--hpo-trials",
        type=int,
        default=50,
        help="Number of Optuna HPO trials per model (default: 50).",
    )
    return parser.parse_args()


# ── Main ───────────────────────────────────────────────────────────────────────

def main() -> None:
    """Entry point."""
    args = parse_args()

    input_dir: Path = args.input_dir
    output_dir: Path = args.output_dir
    results_dir: Path = args.results_dir

    train_path = input_dir / "train_features.csv.gz"
    test_path = input_dir / "test_features.csv.gz"

    for p in (train_path, test_path):
        if not p.exists():
            print(
                f"[ERROR] {p} not found.\n"
                "Run scripts/ml/build_real_training_data.py first.",
                file=sys.stderr,
            )
            sys.exit(1)

    models_dir = output_dir / "models"
    out_results_dir = output_dir / "results"
    models_dir.mkdir(parents=True, exist_ok=True)
    out_results_dir.mkdir(parents=True, exist_ok=True)
    results_dir.mkdir(parents=True, exist_ok=True)

    # ── Load data ──────────────────────────────────────────────────────────────
    print("\n=== Loading data ===", flush=True)
    df_train = load_split(train_path)
    df_test = load_split(test_path)

    X_train, y_train, groups = extract_arrays(df_train)
    X_test, y_test, _ = extract_arrays(df_test)

    n_pos = int(y_train.sum())
    n_neg = int((y_train == 0).sum())
    scale_pos_weight = n_neg / max(n_pos, 1)
    print(f"Train: {len(X_train):,} rows, {n_pos:,} pos, {n_neg:,} neg", flush=True)
    print(f"Test:  {len(X_test):,} rows, {int(y_test.sum()):,} pos", flush=True)
    print(f"Features: {len(FEATURE_NAMES)}", flush=True)

    # ── LightGBM HPO ───────────────────────────────────────────────────────────
    print(f"\n=== LightGBM HPO ({args.hpo_trials} trials) ===", flush=True)
    lgb_study = optuna.create_study(
        direction="maximize",
        sampler=optuna.samplers.TPESampler(seed=42),
    )
    lgb_study.optimize(
        make_lgb_objective(X_train, y_train, groups),
        n_trials=args.hpo_trials,
        show_progress_bar=True,
    )
    lgb_best = lgb_study.best_params
    lgb_best.update({
        "is_unbalance": True,
        "objective": "binary",
        "metric": "average_precision",
        "verbosity": -1,
        "n_jobs": -1,
        "random_state": 42,
    })
    print(f"Best LightGBM AUPRC (HPO CV): {lgb_study.best_value:.4f}", flush=True)

    # ── XGBoost HPO ────────────────────────────────────────────────────────────
    print(f"\n=== XGBoost HPO ({args.hpo_trials} trials) ===", flush=True)
    xgb_study = optuna.create_study(
        direction="maximize",
        sampler=optuna.samplers.TPESampler(seed=42),
    )
    xgb_study.optimize(
        make_xgb_objective(X_train, y_train, groups, scale_pos_weight),
        n_trials=args.hpo_trials,
        show_progress_bar=True,
    )
    xgb_best = xgb_study.best_params
    xgb_best.update({
        "scale_pos_weight": scale_pos_weight,
        "tree_method": "hist",
        "eval_metric": "aucpr",
        "verbosity": 0,
        "n_jobs": -1,
        "random_state": 42,
    })
    print(f"Best XGBoost AUPRC (HPO CV): {xgb_study.best_value:.4f}", flush=True)

    # ── 4-fold CV with best HPO params ─────────────────────────────────────────
    print("\n=== 4-fold GroupKFold CV ===", flush=True)
    all_cv: list[dict] = []

    print("\nLightGBM:", flush=True)
    lgb_cv_model = lgb.LGBMClassifier(**lgb_best)
    all_cv.extend(run_cv(lgb_cv_model, X_train, y_train, groups, "lightgbm"))

    print("\nXGBoost:", flush=True)
    xgb_cv_model = xgb.XGBClassifier(**xgb_best)
    all_cv.extend(run_cv(xgb_cv_model, X_train, y_train, groups, "xgboost"))

    cv_df = pd.DataFrame(all_cv)
    cv_df.to_csv(out_results_dir / "cv_results_v2.csv", index=False)
    print(f"\nCV results written: {out_results_dir / 'cv_results_v2.csv'}", flush=True)

    # ── Final training on full train set ───────────────────────────────────────
    print("\n=== Final training on full training set ===", flush=True)
    lgb_model = lgb.LGBMClassifier(**lgb_best)
    lgb_model.fit(X_train, y_train)
    # Save native LightGBM model so re-export is possible without retraining.
    lgb_model.booster_.save_model(str(out_models_dir / "lightgbm_v2.txt"))

    xgb_model = xgb.XGBClassifier(**xgb_best)
    xgb_model.fit(X_train, y_train)
    xgb_model.save_model(str(out_models_dir / "xgboost_v2.json"))

    # ── Isotonic calibration ───────────────────────────────────────────────────
    print("\n=== Calibrating models (isotonic, 4-fold) ===", flush=True)
    lgb_cal = CalibratedClassifierCV(lgb_model, method="isotonic", cv=4)
    lgb_cal.fit(X_train, y_train)

    xgb_cal = CalibratedClassifierCV(xgb_model, method="isotonic", cv=4)
    xgb_cal.fit(X_train, y_train)

    # ── Test set evaluation ────────────────────────────────────────────────────
    print("\n=== Test set evaluation ===", flush=True)
    prob_lgb = lgb_cal.predict_proba(X_test)[:, 1]
    prob_xgb = xgb_cal.predict_proba(X_test)[:, 1]

    all_test: list[dict] = []
    for model_name, probs in [("lightgbm", prob_lgb), ("xgboost", prob_xgb)]:
        try:
            auroc = float(roc_auc_score(y_test, probs))
            auprc = float(average_precision_score(y_test, probs))
        except ValueError:
            auroc = auprc = float("nan")
        thresh, best_f1 = find_best_f1_threshold(y_test, probs)
        thresh_rows = compute_threshold_metrics(y_test, probs, [0.3, 0.5, 0.7, 0.9])
        print(
            f"  {model_name:<12} AUROC={auroc:.4f}  AUPRC={auprc:.4f}  "
            f"F1={best_f1:.4f}  threshold={thresh:.3f}",
            flush=True,
        )
        entry = {
            "model": model_name,
            "auroc": auroc,
            "auprc": auprc,
            "best_f1": best_f1,
            "best_f1_threshold": thresh,
            "n_pos": int(y_test.sum()),
            "n_neg": int((y_test == 0).sum()),
            "threshold_metrics": thresh_rows,
            "per_vaf_auprc": per_vaf_bin_auprc(df_test, probs, model_name),
        }
        all_test.append(entry)

    with open(out_results_dir / "test_results_v2.json", "w") as fh:
        json.dump(all_test, fh, indent=2)
    print(
        f"\nTest results written: {out_results_dir / 'test_results_v2.json'}",
        flush=True,
    )

    # ── Feature importance ─────────────────────────────────────────────────────
    lgb_imp = pd.DataFrame({
        "feature": FEATURE_NAMES,
        "importance": lgb_model.feature_importances_,
    }).sort_values("importance", ascending=False)
    lgb_imp.to_csv(out_results_dir / "feature_importance_lgb_v2.csv", index=False)

    print("\nTop 20 LightGBM features:", flush=True)
    for _, row in lgb_imp.head(20).iterrows():
        print(f"  {row['feature']:<45s} {row['importance']:.0f}", flush=True)

    # ── Find best F1 threshold on LightGBM for the metadata JSON ──────────────
    lgb_best_thresh, lgb_best_f1 = find_best_f1_threshold(y_test, prob_lgb)

    # ── ONNX export ────────────────────────────────────────────────────────────
    print("\n=== ONNX export ===", flush=True)
    onnx_path = models_dir / "twist-duplex-v2.onnx"
    exported = export_onnx_lgb(lgb_cal, len(FEATURE_NAMES), onnx_path)
    if not exported:
        print("Trying XGBoost ONNX export as fallback.", flush=True)
        export_onnx_xgb(xgb_cal, len(FEATURE_NAMES), onnx_path)

    # ── Model metadata JSON ────────────────────────────────────────────────────
    meta = {
        "version": "twist-duplex-2",
        "feature_names": FEATURE_NAMES,
        "ml_pass_threshold": lgb_best_thresh,
        "variant_class_map": {
            "SNV": 0,
            "Insertion": 1,
            "Deletion": 2,
            "MNV": 3,
            "Complex": 4,
        },
    }
    meta_path = models_dir / "twist-duplex-v2.json"
    with open(meta_path, "w") as fh:
        json.dump(meta, fh, indent=2)
    print(f"Model metadata written: {meta_path}", flush=True)

    # ── Plots ──────────────────────────────────────────────────────────────────
    print("\n=== Saving plots ===", flush=True)
    save_plots(y_test, prob_lgb, results_dir, lgb_imp)

    print("\n=== Training complete ===", flush=True)
    print(
        f"LightGBM test AUPRC: {float(average_precision_score(y_test, prob_lgb)):.4f}  "
        f"ml_pass_threshold: {lgb_best_thresh:.3f}",
        flush=True,
    )


if __name__ == "__main__":
    main()
