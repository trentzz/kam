#!/usr/bin/env python3
"""Generate benchmark figures for the kam paper.

Produces:
  - sensitivity_vs_vaf.pdf   — lollipop chart of sensitivity by VAF and concentration
  - precision_recall.pdf     — scatter plot of precision vs sensitivity
  - molecule_stats.pdf       — molecule count and duplex rate by sample

All figures follow graph-style.md rules:
  - DejaVu Sans, 8pt base font
  - 84mm (single-column) width
  - 600 DPI raster + PDF
  - Thin lines (1.0–1.2pt)
  - Lollipop/dot plots instead of bars
  - No dual y-axes
  - Direct labels where possible
"""

import csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
from pathlib import Path

# ── Paths ─────────────────────────────────────────────────────────────────────
REPO = Path(__file__).resolve().parents[2]
RESULTS = REPO / "benchmarking/results/tables/titration_results.tsv"
OUT_DIR = REPO / "docs/paper/figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ── Style ─────────────────────────────────────────────────────────────────────
MM = 1 / 25.4  # mm to inches
COL_W = 84 * MM  # single column
DPI = 600

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 8,
    "axes.titlesize": 8,
    "axes.labelsize": 8,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 7,
    "lines.linewidth": 1.0,
    "axes.linewidth": 0.6,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.minor.width": 0.4,
    "ytick.minor.width": 0.4,
    "figure.dpi": DPI,
    "savefig.dpi": DPI,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.02,
})

COLOURS = {"5ng": "#2166ac", "15ng": "#4dac26", "30ng": "#d01c8b"}
MARKERS = {"5ng": "o", "15ng": "s", "30ng": "^"}

# ── Load data ─────────────────────────────────────────────────────────────────
def load_results():
    rows = []
    with open(RESULTS) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            rows.append({
                "sample": row["sample"],
                "ng": row["ng"],
                "vaf": float(row["vaf"]),
                "molecules": int(row["molecules"]),
                "duplex": int(row["duplex"]),
                "duplex_pct": float(row["duplex_pct"]),
                "variants_called": int(row["variants_called"]),
                "tp": int(row["tp"]),
                "fp": int(row["fp"]),
                "fn": int(row["fn"]),
                "sensitivity": float(row["sensitivity"]),
                "precision": float(row["precision"]),
                "f1": float(row["f1"]),
                "peak_rss_mb": float(row["peak_rss_mb"]),
                "wall_time_s": float(row["wall_time_s"]),
            })
    return rows


# ── Figure 1: Sensitivity vs VAF (lollipop, one series per concentration) ─────
def plot_sensitivity(rows):
    # Exclude 0% VAF (negative control)
    data = [r for r in rows if r["vaf"] > 0]
    ng_groups = ["5ng", "15ng", "30ng"]
    vafs = sorted({r["vaf"] for r in data})

    fig, ax = plt.subplots(figsize=(COL_W, COL_W * 0.85))

    x_positions = np.arange(len(vafs))
    width = 0.22
    offsets = {"5ng": -width, "15ng": 0, "30ng": width}

    for ng in ng_groups:
        series = {r["vaf"]: r["sensitivity"] for r in data if r["ng"] == ng}
        xpos = x_positions + offsets[ng]
        yvals = [series.get(v, 0.0) for v in vafs]
        ax.vlines(xpos, 0, yvals, color=COLOURS[ng], linewidth=1.0, alpha=0.7)
        ax.plot(xpos, yvals, marker=MARKERS[ng], color=COLOURS[ng],
                markersize=4, linestyle="none", label=ng)

    ax.set_xticks(x_positions)
    ax.set_xticklabels([f"{v:g}" for v in vafs], rotation=30, ha="right")
    ax.set_xlabel("VAF (%)")
    ax.set_ylabel("Sensitivity")
    ax.set_ylim(0, 0.60)
    ax.yaxis.set_major_formatter(ticker.PercentFormatter(xmax=1, decimals=0))
    ax.legend(title="Input", frameon=False, loc="upper left")
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_title("Sensitivity across VAF titration series")

    fig.savefig(OUT_DIR / "sensitivity_vs_vaf.pdf")
    fig.savefig(OUT_DIR / "sensitivity_vs_vaf.png", dpi=DPI)
    plt.close(fig)
    print("Saved sensitivity_vs_vaf")


# ── Figure 2: Precision vs sensitivity scatter ─────────────────────────────────
def plot_precision_recall(rows):
    data = [r for r in rows if r["vaf"] > 0 and r["tp"] > 0]
    ng_groups = ["5ng", "15ng", "30ng"]

    fig, ax = plt.subplots(figsize=(COL_W, COL_W * 0.85))

    for ng in ng_groups:
        series = [r for r in data if r["ng"] == ng]
        xs = [r["sensitivity"] for r in series]
        ys = [r["precision"] for r in series]
        vafs = [r["vaf"] for r in series]
        ax.scatter(xs, ys, marker=MARKERS[ng], color=COLOURS[ng],
                   s=20, label=ng, zorder=3)
        # Label each point with its VAF
        for x, y, v in zip(xs, ys, vafs):
            ax.annotate(f"{v:g}%", (x, y), textcoords="offset points",
                        xytext=(3, 2), fontsize=5.5, color=COLOURS[ng])

    ax.set_xlabel("Sensitivity")
    ax.set_ylabel("Precision")
    ax.xaxis.set_major_formatter(ticker.PercentFormatter(xmax=1, decimals=0))
    ax.yaxis.set_major_formatter(ticker.PercentFormatter(xmax=1, decimals=0))
    ax.set_xlim(0, 0.55)
    ax.set_ylim(0, 1.0)
    ax.legend(title="Input", frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_title("Precision vs sensitivity by concentration and VAF")

    fig.savefig(OUT_DIR / "precision_recall.pdf")
    fig.savefig(OUT_DIR / "precision_recall.png", dpi=DPI)
    plt.close(fig)
    print("Saved precision_recall")


# ── Figure 3: Molecule count and duplex rate ───────────────────────────────────
def plot_molecule_stats(rows):
    # All samples; show by ng group with VAF on x-axis
    ng_groups = ["5ng", "15ng", "30ng"]
    vafs_all = sorted({r["vaf"] for r in rows})

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(COL_W * 1.9, COL_W * 0.75))

    x_positions = np.arange(len(vafs_all))
    width = 0.22
    offsets = {"5ng": -width, "15ng": 0, "30ng": width}

    for ng in ng_groups:
        series = {r["vaf"]: r for r in rows if r["ng"] == ng}
        xpos = x_positions + offsets[ng]
        mols = [series[v]["molecules"] / 1000 if v in series else 0 for v in vafs_all]
        duplex_pct = [series[v]["duplex_pct"] if v in series else 0 for v in vafs_all]
        ax1.vlines(xpos, 0, mols, color=COLOURS[ng], linewidth=1.0, alpha=0.7)
        ax1.plot(xpos, mols, marker=MARKERS[ng], color=COLOURS[ng],
                 markersize=4, linestyle="none", label=ng)
        ax2.vlines(xpos, 0, duplex_pct, color=COLOURS[ng], linewidth=1.0, alpha=0.7)
        ax2.plot(xpos, duplex_pct, marker=MARKERS[ng], color=COLOURS[ng],
                 markersize=4, linestyle="none", label=ng)

    for ax, ylabel, title in [
        (ax1, "Molecules (×1000)", "Assembled molecules per sample"),
        (ax2, "Duplex rate (%)", "Duplex identification rate"),
    ]:
        ax.set_xticks(x_positions)
        ax.set_xticklabels([f"{v:g}" for v in vafs_all], rotation=30, ha="right")
        ax.set_xlabel("VAF (%)")
        ax.set_ylabel(ylabel)
        ax.legend(title="Input", frameon=False, fontsize=6)
        ax.spines[["top", "right"]].set_visible(False)
        ax.set_title(title)

    fig.savefig(OUT_DIR / "molecule_stats.pdf")
    fig.savefig(OUT_DIR / "molecule_stats.png", dpi=DPI)
    plt.close(fig)
    print("Saved molecule_stats")


if __name__ == "__main__":
    rows = load_results()
    plot_sensitivity(rows)
    plot_precision_recall(rows)
    plot_molecule_stats(rows)
    print(f"All figures written to {OUT_DIR}")
