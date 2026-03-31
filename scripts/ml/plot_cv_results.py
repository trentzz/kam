"""Generate ML benchmark figures from CV results and feature importance data.

Outputs five figures to docs/benchmarking/ml/results/figures/:
  1. auroc_comparison.png
  2. auprc_comparison.png
  3. feature_importance.png
  4. pr_curve_summary.png
  5. model_metrics_heatmap.png
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
RESULTS_DIR = REPO_ROOT / "docs/benchmarking/ml/results"
FIGURES_DIR = RESULTS_DIR / "figures"
CV_CSV = RESULTS_DIR / "cv_results.csv"
FI_CSV = RESULTS_DIR / "feature_importance.csv"

FIGURES_DIR.mkdir(parents=True, exist_ok=True)

# Consistent model colours
MODEL_COLOURS = {
    "lightgbm": "#2E86AB",
    "xgboost": "#E84855",
    "baseline_pass": "#6B7280",
}

MODEL_LABELS = {
    "lightgbm": "LightGBM",
    "xgboost": "XGBoost",
    "baseline_pass": "Baseline (PASS)",
}

MODELS = ["lightgbm", "xgboost", "baseline_pass"]


def load_cv(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    df.columns = df.columns.str.strip()
    return df


def load_fi(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    df.columns = df.columns.str.strip()
    return df


def style_axis(ax):
    """Remove top/right spines. Keep it clean."""
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(direction="out", length=4)


def grouped_bar(cv: pd.DataFrame, metric: str, ylabel: str, title: str, fname: str):
    """Grouped bar chart: metric per model across folds with mean ± std error bars."""
    fig, ax = plt.subplots(figsize=(8, 5))

    x = np.arange(5)
    n_models = len(MODELS)
    bar_width = 0.22
    offsets = np.linspace(-(n_models - 1) / 2, (n_models - 1) / 2, n_models) * bar_width

    for i, model in enumerate(MODELS):
        sub = cv[cv["model"] == model].sort_values("fold")
        values = sub[metric].values
        mean_val = values.mean()
        std_val = values.std()

        bars = ax.bar(
            x + offsets[i],
            values,
            width=bar_width,
            color=MODEL_COLOURS[model],
            alpha=0.75,
            label=MODEL_LABELS[model],
            zorder=3,
        )

        # Mean line spanning all fold bars for this model
        span_left = x[0] + offsets[i] - bar_width / 2
        span_right = x[-1] + offsets[i] + bar_width / 2
        ax.hlines(
            mean_val,
            span_left,
            span_right,
            colors=MODEL_COLOURS[model],
            linewidths=2,
            linestyles="--",
            zorder=4,
        )

    ax.set_xlabel("Fold", fontsize=11)
    ax.set_ylabel(ylabel, fontsize=11)
    ax.set_title(title, fontsize=13, fontweight="bold", pad=12)
    ax.set_xticks(x)
    ax.set_xticklabels([f"Fold {i}" for i in range(5)])
    ax.legend(frameon=False, fontsize=10)
    style_axis(ax)
    ax.yaxis.grid(True, linestyle=":", alpha=0.5, zorder=0)
    ax.set_axisbelow(True)

    fig.tight_layout()
    out = FIGURES_DIR / fname
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"Saved {out}")


def plot_feature_importance(fi: pd.DataFrame, fname: str):
    """Horizontal bar chart of top 15 features by LightGBM gain."""
    top = fi.nlargest(15, "importance_gain").sort_values("importance_gain")

    fig, ax = plt.subplots(figsize=(8, 6))

    colours = ["#2E86AB" if v > 0 else "#CBD5E1" for v in top["importance_gain"]]
    ax.barh(top["feature"], top["importance_gain"], color=colours, height=0.6)
    ax.set_xlabel("Gain importance", fontsize=11)
    ax.set_title("LightGBM feature importance (top 15, gain)", fontsize=13, fontweight="bold", pad=12)
    style_axis(ax)
    ax.xaxis.grid(True, linestyle=":", alpha=0.5, zorder=0)
    ax.set_axisbelow(True)

    # Format large numbers with SI suffix
    ax.xaxis.set_major_formatter(
        mticker.FuncFormatter(lambda x, _: f"{x/1e9:.0f}B" if x >= 1e9 else f"{x/1e6:.0f}M")
    )

    fig.tight_layout()
    out = FIGURES_DIR / fname
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"Saved {out}")


def plot_pr_scatter(cv: pd.DataFrame, fname: str):
    """Scatter plot of (recall, precision) per fold per model."""
    fig, ax = plt.subplots(figsize=(7, 5))

    markers = {"lightgbm": "o", "xgboost": "s", "baseline_pass": "^"}

    for model in MODELS:
        sub = cv[cv["model"] == model]
        ax.scatter(
            sub["recall"],
            sub["precision"],
            color=MODEL_COLOURS[model],
            marker=markers[model],
            s=80,
            alpha=0.85,
            label=MODEL_LABELS[model],
            zorder=3,
        )
        # Mean cross-hair
        mr = sub["recall"].mean()
        mp = sub["precision"].mean()
        ax.scatter(
            mr, mp,
            color=MODEL_COLOURS[model],
            marker=markers[model],
            s=200,
            edgecolors="black",
            linewidths=1.2,
            zorder=5,
        )

    ax.set_xlabel("Recall", fontsize=11)
    ax.set_ylabel("Precision", fontsize=11)
    ax.set_title("Precision–recall per fold (large markers = fold mean)", fontsize=13, fontweight="bold", pad=12)
    ax.legend(frameon=False, fontsize=10)
    style_axis(ax)
    ax.yaxis.grid(True, linestyle=":", alpha=0.5, zorder=0)
    ax.xaxis.grid(True, linestyle=":", alpha=0.5, zorder=0)
    ax.set_axisbelow(True)

    fig.tight_layout()
    out = FIGURES_DIR / fname
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"Saved {out}")


def plot_metrics_heatmap(cv: pd.DataFrame, fname: str):
    """Heatmap: rows = models, columns = metrics, values = fold mean."""
    metrics = ["precision", "recall", "f1", "auprc", "auroc"]
    metric_labels = ["Precision", "Recall", "F1", "AUPRC", "AUROC"]

    means = cv.groupby("model")[metrics].mean()
    # Reorder rows
    means = means.loc[MODELS]

    row_labels = [MODEL_LABELS[m] for m in MODELS]

    data = means.values

    fig, ax = plt.subplots(figsize=(8, 3.5))

    # Normalise each column independently for colour mapping
    col_min = data.min(axis=0)
    col_max = data.max(axis=0)
    col_range = np.where(col_max - col_min < 1e-12, 1.0, col_max - col_min)
    norm_data = (data - col_min) / col_range

    im = ax.imshow(norm_data, cmap="Blues", aspect="auto", vmin=0, vmax=1)

    ax.set_xticks(range(len(metrics)))
    ax.set_xticklabels(metric_labels, fontsize=11)
    ax.set_yticks(range(len(MODELS)))
    ax.set_yticklabels(row_labels, fontsize=11)

    # Annotate cells with actual values
    for i in range(len(MODELS)):
        for j in range(len(metrics)):
            val = data[i, j]
            text_colour = "white" if norm_data[i, j] > 0.6 else "black"
            ax.text(j, i, f"{val:.4f}", ha="center", va="center",
                    fontsize=9, color=text_colour)

    ax.set_title("Model metrics (fold mean)", fontsize=13, fontweight="bold", pad=12)
    ax.spines[:].set_visible(False)
    ax.tick_params(length=0)

    fig.tight_layout()
    out = FIGURES_DIR / fname
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"Saved {out}")


def main():
    cv = load_cv(CV_CSV)
    fi = load_fi(FI_CSV)

    grouped_bar(cv, "auroc", "AUROC", "AUROC by model and fold", "auroc_comparison.png")
    grouped_bar(cv, "auprc", "AUPRC", "AUPRC by model and fold", "auprc_comparison.png")
    plot_feature_importance(fi, "feature_importance.png")
    plot_pr_scatter(cv, "pr_curve_summary.png")
    plot_metrics_heatmap(cv, "model_metrics_heatmap.png")

    print("\nAll figures written to", FIGURES_DIR)


if __name__ == "__main__":
    main()
