#!/usr/bin/env python3
"""Evaluate trained Twist duplex models on the held-out test set.

Loads LightGBM, XGBoost, and ensemble models from the models directory and
produces a comprehensive evaluation report: overall metrics, per-VAF-bin
breakdown, per-variant-type breakdown, and curve data for plotting.

Inputs:
  bigdata/experiments/01-ml-twist-duplex/test_features.csv.gz
  bigdata/experiments/01-ml-twist-duplex/models/lightgbm_twist.txt
  bigdata/experiments/01-ml-twist-duplex/models/xgboost_twist.json
  bigdata/experiments/01-ml-twist-duplex/models/ensemble_twist.pkl
  bigdata/experiments/01-ml-twist-duplex/models/xgboost_twist.onnx  (optional)

Outputs (all written to --output-dir):
  evaluation_report.json   — full metrics dict for all models and breakdowns
  roc_data.csv             — fpr, tpr, threshold for each model
  pr_data.csv              — precision, recall, threshold for each model
  confusion_matrix.csv     — counts at the optimal F1 threshold
  vaf_bin_breakdown.csv    — AUPRC, F1, counts per VAF bin per model
  vtype_breakdown.csv      — AUPRC, F1, counts per variant type per model

Usage:
    python3 scripts/ml/evaluate_twist_duplex.py [options]

Run from the repository root.
"""

from __future__ import annotations

import argparse
import json
import logging
import pickle
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.exceptions import UndefinedMetricWarning
from sklearn.metrics import (
    average_precision_score,
    confusion_matrix,
    f1_score,
    precision_recall_curve,
    roc_auc_score,
    roc_curve,
)
from sklearn.preprocessing import LabelEncoder

import lightgbm as lgb
import xgboost as xgb

warnings.filterwarnings("ignore", category=UndefinedMetricWarning)

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

REPO = Path(__file__).resolve().parent.parent.parent

BIGDATA_DIR = REPO / "bigdata" / "experiments" / "01-ml-twist-duplex"
DEFAULT_MODELS_DIR = BIGDATA_DIR / "models"
DEFAULT_TEST_PATH = BIGDATA_DIR / "test_features.csv.gz"

# Small result files tracked in git (mirrors train_twist_duplex.py convention).
DEFAULT_OUTPUT_DIR = REPO / "docs" / "project" / "experiments" / "01-ml-twist-duplex" / "results"

# ---------------------------------------------------------------------------
# Feature definitions (must match train_twist_duplex.py)
# ---------------------------------------------------------------------------

# Columns dropped before building the feature matrix.
DROP_COLS = {
    "target_id", "ref_seq", "alt_seq", "filter",
    "sample_id", "split", "mode", "label",
    "variant_class_true",
}

# Categorical features encoded with LabelEncoder at training time.
CAT_FEATURES = ["variant_class", "sbp_cat", "depth_bucket"]

# VAF bins used for the breakdown report.
VAF_BINS = [
    ("<=0.1%",   0.0,   0.001),
    ("0.1-0.5%", 0.001, 0.005),
    ("0.5-1.0%", 0.005, 0.010),
    ("1.0-5.0%", 0.010, 0.050),
    (">5.0%",    0.050, 1.001),
]

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-8s  %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def load_test_data(path: Path) -> pd.DataFrame:
    """Load the test feature CSV and drop ambiguous rows.

    Removes rows where filter = 'NotTargeted' for consistency with training.

    Args:
        path: Path to the test CSV.gz file.

    Returns:
        Cleaned DataFrame.
    """
    log.info("Loading test data from %s", path)
    df = pd.read_csv(path, low_memory=False)
    log.info("  Loaded %d rows", len(df))

    if "filter" in df.columns:
        before = len(df)
        df = df[df["filter"] != "NotTargeted"]
        removed = before - len(df)
        if removed:
            log.info("  Dropped %d 'NotTargeted' rows", removed)

    n_pos = int((df["label"] == 1).sum())
    n_neg = int((df["label"] == 0).sum())
    log.info("  Positives: %d  Negatives: %d  Rate: %.3f%%",
             n_pos, n_neg, 100 * n_pos / max(len(df), 1))
    return df


def prepare_features(df: pd.DataFrame, feature_cols: list[str] | None = None) -> tuple[pd.DataFrame, list[str]]:
    """Build the feature matrix from a raw DataFrame.

    Encodes categorical columns with LabelEncoder and fills NaN with 0.
    If feature_cols is provided, restricts the matrix to those columns in
    that order (aligns to the model's expected input).

    Args:
        df: Raw DataFrame from load_test_data().
        feature_cols: Optional ordered list of column names from the trained
            model. If None, all non-dropped columns are used.

    Returns:
        Tuple (X, feature_cols) where X is the feature matrix.
    """
    # Encode categoricals in-place before selecting columns.
    df = df.copy()
    for col in CAT_FEATURES:
        if col in df.columns and df[col].dtype == object:
            le = LabelEncoder()
            df[col] = le.fit_transform(df[col].astype(str).fillna("unknown"))

    if feature_cols is not None:
        # Keep only the columns the model expects; fill any missing ones with 0.
        missing = [c for c in feature_cols if c not in df.columns]
        if missing:
            log.warning("  %d feature columns missing from data, filling with 0: %s",
                        len(missing), missing[:5])
        for col in missing:
            df[col] = 0
        X = df[feature_cols].fillna(0)
    else:
        cols = [c for c in df.columns if c not in DROP_COLS]
        X = df[cols].fillna(0)
        feature_cols = cols

    return X, feature_cols


# ---------------------------------------------------------------------------
# Model loading
# ---------------------------------------------------------------------------

def load_lgb_model(models_dir: Path) -> lgb.Booster | None:
    """Load a LightGBM native text model.

    Args:
        models_dir: Directory containing model files.

    Returns:
        Loaded Booster, or None if the file is missing.
    """
    path = models_dir / "lightgbm_twist.txt"
    if not path.exists():
        log.warning("LightGBM model not found: %s", path)
        return None
    log.info("Loading LightGBM model: %s", path)
    return lgb.Booster(model_file=str(path))


def load_xgb_model(models_dir: Path) -> xgb.Booster | None:
    """Load an XGBoost native JSON model.

    Args:
        models_dir: Directory containing model files.

    Returns:
        Loaded Booster, or None if the file is missing.
    """
    path = models_dir / "xgboost_twist.json"
    if not path.exists():
        log.warning("XGBoost model not found: %s", path)
        return None
    log.info("Loading XGBoost model: %s", path)
    booster = xgb.Booster()
    booster.load_model(str(path))
    return booster


def load_ensemble(models_dir: Path) -> dict | None:
    """Load the scikit-learn ensemble pickle.

    The pickle contains a dict with keys 'lgb', 'xgb', and 'feature_cols'.
    Both 'lgb' and 'xgb' are CalibratedClassifierCV wrappers that support
    predict_proba.

    Args:
        models_dir: Directory containing model files.

    Returns:
        Ensemble dict, or None if the file is missing.
    """
    path = models_dir / "ensemble_twist.pkl"
    if not path.exists():
        log.warning("Ensemble pickle not found: %s", path)
        return None
    log.info("Loading ensemble pickle: %s", path)
    with open(path, "rb") as fh:
        return pickle.load(fh)


def verify_onnx(models_dir: Path, X: np.ndarray, prob_ref: np.ndarray) -> None:
    """Verify the ONNX export against the reference probabilities.

    Compares ONNX output to prob_ref (XGBoost sklearn probabilities) and
    logs the mean absolute difference. Skips gracefully if onnxruntime is not
    available or the model file is missing.

    Args:
        models_dir: Directory containing model files.
        X: Feature matrix (float32 input expected by ONNX).
        prob_ref: Reference probability array to compare against.
    """
    path = models_dir / "xgboost_twist.onnx"
    if not path.exists():
        log.info("ONNX model not found, skipping verification.")
        return

    try:
        import onnxruntime as ort
    except ImportError:
        log.warning("onnxruntime not available — skipping ONNX verification.")
        return

    log.info("Verifying ONNX export: %s", path)
    sess = ort.InferenceSession(str(path))
    input_name = sess.get_inputs()[0].name
    # ONNX models expect float32.
    out = sess.run(None, {input_name: X.astype(np.float32)})
    # Output is typically [labels, probabilities_dict]; extract class-1 prob.
    if isinstance(out[1], list):
        prob_onnx = np.array([d[1] for d in out[1]], dtype=float)
    else:
        prob_onnx = np.array(out[1])[:, 1]
    mae = float(np.mean(np.abs(prob_onnx - prob_ref)))
    log.info("  ONNX vs reference mean absolute error: %.6f", mae)


# ---------------------------------------------------------------------------
# Metrics
# ---------------------------------------------------------------------------

def best_f1_threshold(y_true: np.ndarray, y_prob: np.ndarray) -> tuple[float, float]:
    """Find the threshold that maximises F1 on the precision-recall curve.

    Args:
        y_true: Ground-truth binary labels.
        y_prob: Predicted probabilities for the positive class.

    Returns:
        Tuple (best_f1, threshold).
    """
    prec, rec, thresholds = precision_recall_curve(y_true, y_prob)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        f1_scores = 2 * prec * rec / (prec + rec + 1e-9)
    idx = int(np.argmax(f1_scores))
    # precision_recall_curve returns one more value than thresholds.
    threshold = float(thresholds[min(idx, len(thresholds) - 1)])
    return float(f1_scores[idx]), threshold


def specificity_at_sensitivity(
    fpr_arr: np.ndarray,
    tpr_arr: np.ndarray,
    target_sensitivity: float,
) -> float:
    """Return specificity (1 - FPR) at the given sensitivity level.

    Interpolates along the ROC curve.

    Args:
        fpr_arr: False positive rates from roc_curve().
        tpr_arr: True positive rates from roc_curve().
        target_sensitivity: The TPR level to anchor at (e.g. 0.90 or 0.95).

    Returns:
        Specificity value in [0, 1].
    """
    # Interpolate FPR as a function of TPR.
    fpr_at_target = float(np.interp(target_sensitivity, tpr_arr, fpr_arr))
    return 1.0 - fpr_at_target


def compute_overall_metrics(
    y_true: np.ndarray,
    y_prob: np.ndarray,
    model_name: str,
    threshold: float | None = None,
) -> dict:
    """Compute the full set of overall evaluation metrics.

    Metrics include AUROC, AUPRC, best F1 and its threshold, sensitivity at
    1% and 5% FPR, and specificity at 90% and 95% sensitivity.

    Args:
        y_true: Ground-truth binary labels.
        y_prob: Predicted probabilities for the positive class.
        model_name: Label to attach to the result dict.
        threshold: Decision threshold. If None, uses the best-F1 threshold.

    Returns:
        Dict of metric values.
    """
    try:
        auroc = float(roc_auc_score(y_true, y_prob))
        auprc = float(average_precision_score(y_true, y_prob))
    except ValueError:
        auroc = auprc = float("nan")

    fpr_arr, tpr_arr, _ = roc_curve(y_true, y_prob)
    sens_at_1pct_fpr = float(np.interp(0.01, fpr_arr, tpr_arr))
    sens_at_5pct_fpr = float(np.interp(0.05, fpr_arr, tpr_arr))
    spec_at_90pct_sens = specificity_at_sensitivity(fpr_arr, tpr_arr, 0.90)
    spec_at_95pct_sens = specificity_at_sensitivity(fpr_arr, tpr_arr, 0.95)

    f1_best, f1_threshold = best_f1_threshold(y_true, y_prob)
    decision_threshold = threshold if threshold is not None else f1_threshold

    y_pred = (y_prob >= decision_threshold).astype(int)
    f1_at_threshold = float(f1_score(y_true, y_pred, zero_division=0))

    return {
        "model": model_name,
        "auroc": auroc,
        "auprc": auprc,
        "best_f1": f1_best,
        "best_f1_threshold": f1_threshold,
        "decision_threshold": decision_threshold,
        "f1_at_threshold": f1_at_threshold,
        "sens_at_1pct_fpr": sens_at_1pct_fpr,
        "sens_at_5pct_fpr": sens_at_5pct_fpr,
        "spec_at_90pct_sens": spec_at_90pct_sens,
        "spec_at_95pct_sens": spec_at_95pct_sens,
        "n_pos": int(y_true.sum()),
        "n_neg": int((y_true == 0).sum()),
    }


def compute_vaf_bin_metrics(
    y_true: np.ndarray,
    y_prob: np.ndarray,
    vaf: np.ndarray,
    model_name: str,
) -> list[dict]:
    """Compute AUPRC and best F1 for each VAF bin.

    Bins with fewer than 10 samples are skipped.

    Args:
        y_true: Ground-truth binary labels.
        y_prob: Predicted probabilities.
        vaf: Per-row VAF values (from vaf_target column).
        model_name: Model identifier.

    Returns:
        List of metric dicts, one per populated VAF bin.
    """
    results = []
    for bin_name, vaf_low, vaf_high in VAF_BINS:
        mask = (vaf >= vaf_low) & (vaf < vaf_high)
        n = int(mask.sum())
        n_pos = int(y_true[mask].sum())
        n_neg = n - n_pos
        if n < 10:
            continue

        try:
            auprc = float(average_precision_score(y_true[mask], y_prob[mask]))
        except ValueError:
            auprc = float("nan")

        f1_best, _ = best_f1_threshold(y_true[mask], y_prob[mask])

        results.append({
            "model": model_name,
            "vaf_bin": bin_name,
            "n_positive": n_pos,
            "n_negative": n_neg,
            "auprc": auprc,
            "best_f1": f1_best,
        })
    return results


def compute_vtype_metrics(
    y_true: np.ndarray,
    y_prob: np.ndarray,
    vtypes: np.ndarray,
    model_name: str,
) -> list[dict]:
    """Compute AUPRC and best F1 for each variant type.

    Types with fewer than 10 samples are skipped.

    Args:
        y_true: Ground-truth binary labels.
        y_prob: Predicted probabilities.
        vtypes: Per-row variant type strings (from variant_class_true column).
        model_name: Model identifier.

    Returns:
        List of metric dicts, one per populated variant type.
    """
    results = []
    for vtype in sorted(set(vtypes)):
        mask = vtypes == vtype
        n = int(mask.sum())
        n_pos = int(y_true[mask].sum())
        n_neg = n - n_pos
        if n < 10:
            continue

        try:
            auprc = float(average_precision_score(y_true[mask], y_prob[mask]))
        except ValueError:
            auprc = float("nan")

        f1_best, _ = best_f1_threshold(y_true[mask], y_prob[mask])

        results.append({
            "model": model_name,
            "vtype": vtype,
            "n_positive": n_pos,
            "n_negative": n_neg,
            "auprc": auprc,
            "best_f1": f1_best,
        })
    return results


# ---------------------------------------------------------------------------
# Curve data
# ---------------------------------------------------------------------------

def build_roc_rows(y_true: np.ndarray, y_prob: np.ndarray, model_name: str) -> list[dict]:
    """Build rows for the ROC curve CSV.

    Args:
        y_true: Ground-truth labels.
        y_prob: Predicted probabilities.
        model_name: Model identifier.

    Returns:
        List of dicts with keys model, fpr, tpr, threshold.
    """
    fpr, tpr, thresholds = roc_curve(y_true, y_prob)
    rows = []
    for f, t, th in zip(fpr, tpr, thresholds):
        rows.append({"model": model_name, "fpr": float(f), "tpr": float(t), "threshold": float(th)})
    return rows


def build_pr_rows(y_true: np.ndarray, y_prob: np.ndarray, model_name: str) -> list[dict]:
    """Build rows for the precision-recall curve CSV.

    Args:
        y_true: Ground-truth labels.
        y_prob: Predicted probabilities.
        model_name: Model identifier.

    Returns:
        List of dicts with keys model, precision, recall, threshold.
        The last entry has no threshold (precision_recall_curve convention).
    """
    prec, rec, thresholds = precision_recall_curve(y_true, y_prob)
    rows = []
    for i, (p, r) in enumerate(zip(prec, rec)):
        th = float(thresholds[i]) if i < len(thresholds) else float("nan")
        rows.append({"model": model_name, "precision": float(p), "recall": float(r), "threshold": th})
    return rows


def build_confusion_rows(
    y_true: np.ndarray,
    y_prob: np.ndarray,
    threshold: float,
    model_name: str,
) -> list[dict]:
    """Build rows for the confusion matrix CSV at a given threshold.

    Args:
        y_true: Ground-truth labels.
        y_prob: Predicted probabilities.
        threshold: Decision threshold.
        model_name: Model identifier.

    Returns:
        List of dicts with keys model, threshold, true_label, pred_label, count.
    """
    y_pred = (y_prob >= threshold).astype(int)
    cm = confusion_matrix(y_true, y_pred)
    rows = []
    for i, true_label in enumerate([0, 1]):
        for j, pred_label in enumerate([0, 1]):
            rows.append({
                "model": model_name,
                "threshold": threshold,
                "true_label": true_label,
                "pred_label": pred_label,
                "count": int(cm[i, j]),
            })
    return rows


# ---------------------------------------------------------------------------
# Prediction helpers
# ---------------------------------------------------------------------------

def predict_lgb(model: lgb.Booster, X: pd.DataFrame) -> np.ndarray:
    """Run LightGBM native Booster prediction.

    The native Booster.predict() returns P(class 1) directly.

    Args:
        model: Loaded lgb.Booster.
        X: Feature matrix as a DataFrame.

    Returns:
        1-D probability array.
    """
    return model.predict(X.values)


def predict_xgb(model: xgb.Booster, X: pd.DataFrame) -> np.ndarray:
    """Run XGBoost native Booster prediction via DMatrix.

    Args:
        model: Loaded xgb.Booster.
        X: Feature matrix as a DataFrame.

    Returns:
        1-D probability array.
    """
    dmat = xgb.DMatrix(X.values, feature_names=list(X.columns))
    return model.predict(dmat)


def predict_ensemble(ensemble: dict, X: pd.DataFrame) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Run ensemble prediction, returning per-model and averaged probabilities.

    The ensemble pickle stores calibrated LGB and XGB classifiers under keys
    'lgb' and 'xgb', plus the expected 'feature_cols' list.

    Args:
        ensemble: Dict loaded from ensemble_twist.pkl.
        X: Feature matrix as a DataFrame aligned to feature_cols.

    Returns:
        Tuple (prob_lgb, prob_xgb, prob_ensemble).
    """
    prob_lgb = ensemble["lgb"].predict_proba(X.values)[:, 1]
    prob_xgb = ensemble["xgb"].predict_proba(X.values)[:, 1]
    prob_ensemble = 0.5 * prob_lgb + 0.5 * prob_xgb
    return prob_lgb, prob_xgb, prob_ensemble


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Evaluate trained Twist duplex models on the held-out test set.",
    )
    parser.add_argument(
        "--models-dir",
        type=Path,
        default=DEFAULT_MODELS_DIR,
        help=f"Directory containing trained model files (default: {DEFAULT_MODELS_DIR}).",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help=f"Directory to write evaluation outputs (default: {DEFAULT_OUTPUT_DIR}).",
    )
    parser.add_argument(
        "--test-features",
        type=Path,
        default=DEFAULT_TEST_PATH,
        help=f"Path to test_features.csv.gz (default: {DEFAULT_TEST_PATH}).",
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=None,
        help="Decision threshold. If omitted, the best-F1 threshold is used per model.",
    )
    return parser.parse_args()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    """Entry point."""
    args = parse_args()

    if not args.test_features.exists():
        log.error(
            "%s not found.\nRun scripts/ml/build_twist_duplex_features.py first.",
            args.test_features,
        )
        sys.exit(1)

    args.output_dir.mkdir(parents=True, exist_ok=True)

    # ── Load models ───────────────────────────────────────────────────────────

    log.info("=== Loading models ===")
    lgb_model = load_lgb_model(args.models_dir)
    xgb_model = load_xgb_model(args.models_dir)
    ensemble = load_ensemble(args.models_dir)

    if lgb_model is None and xgb_model is None and ensemble is None:
        log.error("No models could be loaded from %s. Nothing to evaluate.", args.models_dir)
        sys.exit(1)

    # Determine expected feature columns from the ensemble pickle if available,
    # so that all models use a consistent feature set.
    ensemble_feature_cols: list[str] | None = None
    if ensemble is not None:
        ensemble_feature_cols = ensemble.get("feature_cols")

    # ── Load and prepare data ─────────────────────────────────────────────────

    log.info("=== Loading test data ===")
    df_test = load_test_data(args.test_features)
    y_true = df_test["label"].values
    vaf = df_test["vaf_target"].fillna(0).values if "vaf_target" in df_test.columns else np.zeros(len(df_test))
    vtypes = df_test["variant_class_true"].fillna("unknown").values if "variant_class_true" in df_test.columns else np.full(len(df_test), "unknown")

    # Prepare feature matrices for each model.  The ensemble feature list takes
    # priority; native boosters fall back to that list if available.
    X_df, inferred_cols = prepare_features(df_test, feature_cols=ensemble_feature_cols)
    log.info("Feature matrix: %d rows x %d columns", len(X_df), len(inferred_cols))

    # ── Build predictions ─────────────────────────────────────────────────────

    log.info("=== Running predictions ===")
    probs: dict[str, np.ndarray] = {}

    # Ensemble (calibrated sklearn wrappers) — most reliable probabilities.
    if ensemble is not None:
        prob_lgb_cal, prob_xgb_cal, prob_ens = predict_ensemble(ensemble, X_df)
        probs["lgb_calibrated"] = prob_lgb_cal
        probs["xgb_calibrated"] = prob_xgb_cal
        probs["ensemble"] = prob_ens
        log.info("  Ensemble predictions done.")

        # Use calibrated XGBoost output as the ONNX reference.
        verify_onnx(args.models_dir, X_df.values, prob_xgb_cal)

    # Native LightGBM Booster (uncalibrated, but directly interpretable).
    if lgb_model is not None:
        probs["lightgbm_native"] = predict_lgb(lgb_model, X_df)
        log.info("  LightGBM native predictions done.")

    # Native XGBoost Booster (uncalibrated).
    if xgb_model is not None:
        probs["xgboost_native"] = predict_xgb(xgb_model, X_df)
        log.info("  XGBoost native predictions done.")

    # ── Compute metrics ───────────────────────────────────────────────────────

    log.info("=== Computing metrics ===")
    overall_metrics: list[dict] = []
    vaf_bin_rows: list[dict] = []
    vtype_rows: list[dict] = []
    roc_rows: list[dict] = []
    pr_rows: list[dict] = []
    confusion_rows: list[dict] = []

    for model_name, y_prob in probs.items():
        log.info("  %s", model_name)
        m = compute_overall_metrics(y_true, y_prob, model_name, threshold=args.threshold)
        overall_metrics.append(m)

        threshold = m["decision_threshold"]
        log.info(
            "    AUROC=%.4f  AUPRC=%.4f  F1=%.4f  threshold=%.4f",
            m["auroc"], m["auprc"], m["best_f1"], threshold,
        )

        vaf_bin_rows.extend(compute_vaf_bin_metrics(y_true, y_prob, vaf, model_name))
        vtype_rows.extend(compute_vtype_metrics(y_true, y_prob, vtypes, model_name))
        roc_rows.extend(build_roc_rows(y_true, y_prob, model_name))
        pr_rows.extend(build_pr_rows(y_true, y_prob, model_name))
        confusion_rows.extend(build_confusion_rows(y_true, y_prob, threshold, model_name))

    # ── Write outputs ─────────────────────────────────────────────────────────

    log.info("=== Writing outputs to %s ===", args.output_dir)

    # Full metrics report.
    report = {
        "overall": overall_metrics,
        "vaf_bin_breakdown": vaf_bin_rows,
        "vtype_breakdown": vtype_rows,
    }
    report_path = args.output_dir / "evaluation_report.json"
    with open(report_path, "w") as fh:
        json.dump(report, fh, indent=2)
    log.info("  Written: %s", report_path)

    # Curve data CSVs.
    roc_path = args.output_dir / "roc_data.csv"
    pd.DataFrame(roc_rows).to_csv(roc_path, index=False)
    log.info("  Written: %s", roc_path)

    pr_path = args.output_dir / "pr_data.csv"
    pd.DataFrame(pr_rows).to_csv(pr_path, index=False)
    log.info("  Written: %s", pr_path)

    confusion_path = args.output_dir / "confusion_matrix.csv"
    pd.DataFrame(confusion_rows).to_csv(confusion_path, index=False)
    log.info("  Written: %s", confusion_path)

    vaf_path = args.output_dir / "vaf_bin_breakdown.csv"
    pd.DataFrame(vaf_bin_rows).to_csv(vaf_path, index=False)
    log.info("  Written: %s", vaf_path)

    vtype_path = args.output_dir / "vtype_breakdown.csv"
    pd.DataFrame(vtype_rows).to_csv(vtype_path, index=False)
    log.info("  Written: %s", vtype_path)

    log.info("=== Evaluation complete ===")


if __name__ == "__main__":
    main()
