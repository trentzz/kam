#!/usr/bin/env python3
"""Generate benchmark figures for the kam paper.

Single-run mode (default):
  python3 plot_results.py [--results PATH]

  Produces:
    sensitivity_vs_vaf.pdf   — grouped bar chart of sensitivity by VAF and concentration
    precision_recall.pdf     — scatter plot of precision vs sensitivity
    molecule_stats.pdf       — molecule count and duplex rate by sample
    confusion_matrix.pdf     — TP/FP/FN counts at 2% VAF by concentration
    fn_analysis.pdf          — missed variants by type across VAF titration

Compare mode:
  python3 plot_results.py --compare 250k:titration_results_250kreads.tsv \\
                                    500k:titration_results_500kreads.tsv \\
                                    1m:titration_results_1mreads.tsv \\
                                    2m:titration_results_2mreads.tsv

  Produces (averaged across all three concentrations):
    compare_sensitivity_vs_vaf.pdf        — overall sensitivity by read depth
    compare_sensitivity_snv_vs_vaf.pdf    — SNV sensitivity by read depth
    compare_sensitivity_indel_vs_vaf.pdf  — indel sensitivity by read depth
    compare_precision_vs_vaf.pdf          — precision by read depth

All figures follow graph-style.md rules:
  - DejaVu Sans, 8pt base font
  - 84mm (single-column) width
  - 600 DPI raster + PDF
  - Thin lines (1.0–1.2pt)
  - Grouped vertical bar charts
  - No dual y-axes
  - Direct labels where possible
"""

import argparse
import csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
from pathlib import Path

# ── Paths ─────────────────────────────────────────────────────────────────────
REPO = Path(__file__).resolve().parents[2]
_DEFAULT_RESULTS = REPO / "benchmarking/results/tables_v6a/titration_results_2mreads.tsv"
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

# Colours for read-depth comparison plots (colour-blind friendly).
DEPTH_COLOURS = ["#1b7837", "#762a83", "#e08214", "#2166ac"]

# ── Load data ─────────────────────────────────────────────────────────────────
def load_results(path):
    """Load a titration results TSV. Skips rows with missing/empty metric fields."""
    rows = []
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            # Skip failed/skipped rows (exit_code == -1 or empty metrics).
            if not row.get("molecules"):
                continue
            try:
                rows.append({
                    "sample": row["sample"],
                    "ng": row["ng"],
                    "vaf": float(row["vaf"]),
                    "molecules": int(row["molecules"]),
                    "duplex": int(row["duplex"]),
                    "duplex_pct": float(row["duplex_pct"]),
                    "variants_called": int(row["variants_called"]),
                    # overall
                    "tp": int(row["tp"]),
                    "fp": int(row["fp"]),
                    "fn": int(row["fn"]),
                    "sensitivity": float(row["sensitivity"]),
                    "precision": float(row["precision"]),
                    "f1": float(row["f1"]),
                    # SNV
                    "snv_tp": int(row["snv_tp"]),
                    "snv_fp": int(row["snv_fp"]),
                    "snv_fn": int(row["snv_fn"]),
                    "snv_sensitivity": float(row["snv_sensitivity"]),
                    "snv_precision": float(row["snv_precision"]),
                    # indel
                    "indel_tp": int(row["indel_tp"]),
                    "indel_fp": int(row["indel_fp"]),
                    "indel_fn": int(row["indel_fn"]),
                    "indel_sensitivity": float(row["indel_sensitivity"]),
                    "indel_precision": float(row["indel_precision"]),
                    # runtime summary
                    "peak_rss_mb": float(row["peak_rss_mb"]),
                    "wall_time_s": float(row["wall_time_s"]),
                    # per-stage timing (ms) and RSS (MB)
                    "t_assemble_ms":  float(row.get("t_assemble_ms",  0) or 0),
                    "t_index_ms":     float(row.get("t_index_ms",     0) or 0),
                    "t_pathfind_ms":  float(row.get("t_pathfind_ms",  0) or 0),
                    "t_call_ms":      float(row.get("t_call_ms",      0) or 0),
                    "t_output_ms":    float(row.get("t_output_ms",    0) or 0),
                    "rss_assemble_mb": float(row.get("rss_assemble_mb", 0) or 0),
                    "rss_index_mb":    float(row.get("rss_index_mb",    0) or 0),
                    "rss_pathfind_mb": float(row.get("rss_pathfind_mb", 0) or 0),
                    "rss_call_mb":     float(row.get("rss_call_mb",     0) or 0),
                    "rss_output_mb":   float(row.get("rss_output_mb",   0) or 0),
                })
            except (ValueError, KeyError):
                continue
    return rows


# ── Single-run figures ─────────────────────────────────────────────────────────

def plot_sensitivity(rows):
    """Grouped vertical bar chart of overall sensitivity by VAF and concentration."""
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
        ax.bar(xpos, yvals, width=width, color=COLOURS[ng], label=ng,
               linewidth=0.5, edgecolor="white")

    ax.set_xticks(x_positions)
    ax.set_xticklabels([f"{v:g}" for v in vafs], rotation=30, ha="right")
    ax.set_xlabel("VAF (%)")
    ax.set_ylabel("Sensitivity")
    ax.set_ylim(0, 0.70)
    ax.yaxis.set_major_formatter(ticker.PercentFormatter(xmax=1, decimals=0))
    ax.legend(title="Input", frameon=False, loc="upper left")
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_title("Sensitivity across VAF titration series")

    fig.savefig(OUT_DIR / "sensitivity_vs_vaf.pdf")
    fig.savefig(OUT_DIR / "sensitivity_vs_vaf.png", dpi=DPI)
    plt.close(fig)
    print("Saved sensitivity_vs_vaf")


def plot_precision_recall(rows):
    """Scatter plot of precision vs sensitivity, labelled by VAF."""
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


def plot_molecule_stats(rows):
    """Two-panel grouped bar chart: molecule count and duplex rate by VAF and concentration."""
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
        ax1.bar(xpos, mols, width=width, color=COLOURS[ng], label=ng,
                linewidth=0.5, edgecolor="white")
        ax2.bar(xpos, duplex_pct, width=width, color=COLOURS[ng], label=ng,
                linewidth=0.5, edgecolor="white")

    for ax, ylabel, title in [
        (ax1, "Molecules (×1000)", "Assembled molecules per sample"),
        (ax2, "Duplex rate (%)", "Duplex identification rate"),
    ]:
        ax.set_xticks(x_positions)
        ax.set_xticklabels([f"{v:g}" for v in vafs_all], rotation=30, ha="right")
        ax.set_xlabel("VAF (%)")
        ax.set_ylabel(ylabel)
        ax.set_ylim(bottom=0)
        ax.legend(title="Input", frameon=False, fontsize=6)
        ax.spines[["top", "right"]].set_visible(False)
        ax.set_title(title)

    fig.savefig(OUT_DIR / "molecule_stats.pdf")
    fig.savefig(OUT_DIR / "molecule_stats.png", dpi=DPI)
    plt.close(fig)
    print("Saved molecule_stats")


# ── Figure 4: SNV vs indel sensitivity ────────────────────────────────────────

def plot_sensitivity_by_type(rows):
    """Two-panel grouped bar chart: SNV sensitivity (left) and indel sensitivity (right)."""
    data = [r for r in rows if r["vaf"] > 0]
    ng_groups = ["5ng", "15ng", "30ng"]
    vafs = sorted({r["vaf"] for r in data})

    x_positions = np.arange(len(vafs))
    width = 0.22
    offsets = {"5ng": -width, "15ng": 0, "30ng": width}

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(COL_W * 1.9, COL_W * 0.85),
                                   sharey=False)

    for ax, metric, title in [
        (ax1, "snv_sensitivity",   "SNV sensitivity"),
        (ax2, "indel_sensitivity", "Indel sensitivity"),
    ]:
        for ng in ng_groups:
            series = {r["vaf"]: r[metric] for r in data if r["ng"] == ng}
            xpos = x_positions + offsets[ng]
            yvals = [series.get(v, 0.0) for v in vafs]
            ax.bar(xpos, yvals, width=width, color=COLOURS[ng], label=ng,
                   linewidth=0.5, edgecolor="white")
        ax.set_xticks(x_positions)
        ax.set_xticklabels([f"{v:g}" for v in vafs], rotation=30, ha="right")
        ax.set_xlabel("VAF (%)")
        ax.set_ylabel("Sensitivity")
        ax.set_ylim(0, 1.0)
        ax.yaxis.set_major_formatter(ticker.PercentFormatter(xmax=1, decimals=0))
        ax.legend(title="Input", frameon=False, loc="upper left", fontsize=6)
        ax.spines[["top", "right"]].set_visible(False)
        ax.set_title(title)

    fig.savefig(OUT_DIR / "sensitivity_by_type.pdf")
    fig.savefig(OUT_DIR / "sensitivity_by_type.png", dpi=DPI)
    plt.close(fig)
    print("Saved sensitivity_by_type")


# ── Figure 5: False positives by VAF ──────────────────────────────────────────

def plot_fp_by_vaf(rows):
    """Grouped vertical bar chart of false positive count by VAF and concentration."""
    data = [r for r in rows if r["vaf"] > 0]
    ng_groups = ["5ng", "15ng", "30ng"]
    vafs = sorted({r["vaf"] for r in data})

    x_positions = np.arange(len(vafs))
    width = 0.22
    offsets = {"5ng": -width, "15ng": 0, "30ng": width}

    fig, ax = plt.subplots(figsize=(COL_W, COL_W * 0.85))

    for ng in ng_groups:
        series = {r["vaf"]: r["fp"] for r in data if r["ng"] == ng}
        xpos = x_positions + offsets[ng]
        yvals = [series.get(v, 0) for v in vafs]
        ax.bar(xpos, yvals, width=width, color=COLOURS[ng], label=ng,
               linewidth=0.5, edgecolor="white")

    ax.set_xticks(x_positions)
    ax.set_xticklabels([f"{v:g}" for v in vafs], rotation=30, ha="right")
    ax.set_xlabel("VAF (%)")
    ax.set_ylabel("False positives")
    ax.set_ylim(bottom=0)
    ax.legend(title="Input", frameon=False, loc="upper right")
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_title("False positive calls by VAF and concentration")

    fig.savefig(OUT_DIR / "fp_by_vaf.pdf")
    fig.savefig(OUT_DIR / "fp_by_vaf.png", dpi=DPI)
    plt.close(fig)
    print("Saved fp_by_vaf")


# ── Figure 6: Per-stage runtime and memory breakdown ──────────────────────────

def plot_stage_breakdown(rows):
    """Two-panel vertical bar chart: per-stage timing and peak RSS.

    Each bar shows the median value for that stage across all samples.
    """
    stages = ["assemble", "index", "pathfind", "call", "output"]
    stage_labels = ["Assemble", "Index", "Pathfind", "Call", "Output"]

    timing_key = lambda s: f"t_{s}_ms"
    rss_key    = lambda s: f"rss_{s}_mb"

    # Collect all values per stage.
    timing_vals = {s: [r[timing_key(s)] for r in rows
                       if timing_key(s) in r and r[timing_key(s)] > 0]
                   for s in stages}
    rss_vals    = {s: [r[rss_key(s)] for r in rows
                       if rss_key(s) in r and r[rss_key(s)] > 0]
                   for s in stages}

    # Check we have per-stage data.
    if not any(timing_vals[s] for s in stages):
        print("Skipping stage_breakdown: no per-stage timing data")
        return

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(COL_W * 1.9, COL_W * 0.75))

    x_positions = np.arange(len(stages))

    for ax, vals_dict, ylabel, title in [
        (ax1, timing_vals, "Wall time (ms)", "Per-stage runtime"),
        (ax2, rss_vals,    "Peak RSS (MB)",  "Per-stage memory"),
    ]:
        medians = []
        for s in stages:
            vals = vals_dict[s]
            medians.append(np.median(vals) if vals else 0.0)

        ax.bar(x_positions, medians, width=0.55, color="#2166ac",
               linewidth=0.5, edgecolor="white")

        ax.set_xticks(x_positions)
        ax.set_xticklabels(stage_labels, rotation=30, ha="right")
        ax.set_ylabel(ylabel)
        ax.set_ylim(bottom=0)
        ax.spines[["top", "right"]].set_visible(False)
        ax.set_title(title)

    fig.savefig(OUT_DIR / "stage_breakdown.pdf")
    fig.savefig(OUT_DIR / "stage_breakdown.png", dpi=DPI)
    plt.close(fig)
    print("Saved stage_breakdown")


# ── Figure 7: Confusion matrix at 2% VAF ──────────────────────────────────────

def plot_confusion_matrix(rows):
    """Grouped bar chart of TP, FP, and FN counts at 2% VAF by concentration.

    Shows three bars per concentration group: TP (green), FP (pink), FN (red-orange).
    A dashed horizontal line marks the total truth variant count (375).
    """
    # Filter to 2% VAF only.
    target_vaf = 2.0
    data = [r for r in rows if abs(r["vaf"] - target_vaf) < 1e-6]

    if not data:
        print("Skipping confusion_matrix: no rows at 2% VAF")
        return

    ng_groups = ["5ng", "15ng", "30ng"]
    total_truth = 375

    # Colours for TP, FP, FN.
    tp_colour = "#1b7837"
    fp_colour = "#d01c8b"
    fn_colour = "#d6604d"

    fig, ax = plt.subplots(figsize=(COL_W, COL_W * 0.85))

    # Three concentrations on x-axis; within each, three bars side by side.
    x_positions = np.arange(len(ng_groups))
    bar_width = 0.22
    offsets = {"tp": -bar_width, "fp": 0, "fn": bar_width}
    metric_colours = {"tp": tp_colour, "fp": fp_colour, "fn": fn_colour}
    metric_labels  = {"tp": "TP", "fp": "FP", "fn": "FN"}

    # Build lookup: ng -> row.
    by_ng = {r["ng"]: r for r in data}

    for metric in ("tp", "fp", "fn"):
        yvals = [by_ng[ng][metric] if ng in by_ng else 0 for ng in ng_groups]
        xpos = x_positions + offsets[metric]
        bars = ax.bar(xpos, yvals, width=bar_width,
                      color=metric_colours[metric], label=metric_labels[metric],
                      linewidth=0.5, edgecolor="white")
        # Value labels on top of each bar.
        for bar, val in zip(bars, yvals):
            if val > 0:
                ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 2,
                        str(val), ha="center", va="bottom", fontsize=5.5)

    # Dashed line for total truth variant count.
    ax.axhline(total_truth, color="#555555", linewidth=0.8, linestyle="--",
               label=f"Truth total ({total_truth})")

    ax.set_xticks(x_positions)
    ax.set_xticklabels(ng_groups)
    ax.set_xlabel("Input concentration")
    ax.set_ylabel("Variant count")
    ax.set_ylim(bottom=0)
    ax.legend(frameon=False, loc="upper right", fontsize=6)
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_title("True/false positive and false negative counts at 2% VAF")

    fig.savefig(OUT_DIR / "confusion_matrix.pdf")
    fig.savefig(OUT_DIR / "confusion_matrix.png", dpi=DPI)
    plt.close(fig)
    print("Saved confusion_matrix")


# ── Figure 8: FN breakdown by variant type ────────────────────────────────────

def plot_fn_analysis(rows):
    """Grouped bar chart of missed variants (FN) by type across VAF levels.

    SNV FN and indel FN values are averaged across concentrations at each VAF level.
    Total truth: 205 SNVs, 170 indels.
    """
    from collections import defaultdict

    data = [r for r in rows if r["vaf"] > 0]
    vafs = sorted({r["vaf"] for r in data})

    # Average SNV FN and indel FN across concentrations at each VAF.
    snv_fn_by_vaf   = defaultdict(list)
    indel_fn_by_vaf = defaultdict(list)
    for r in data:
        snv_fn_by_vaf[r["vaf"]].append(r["snv_fn"])
        indel_fn_by_vaf[r["vaf"]].append(r["indel_fn"])

    snv_fn_mean   = [np.mean(snv_fn_by_vaf[v])   for v in vafs]
    indel_fn_mean = [np.mean(indel_fn_by_vaf[v]) for v in vafs]

    snv_colour   = "#2166ac"
    indel_colour = "#d01c8b"

    fig, ax = plt.subplots(figsize=(COL_W, COL_W * 0.85))

    x_positions = np.arange(len(vafs))
    bar_width = 0.30
    offsets = {"snv": -bar_width / 2, "indel": bar_width / 2}

    ax.bar(x_positions + offsets["snv"],   snv_fn_mean,   width=bar_width,
           color=snv_colour,   label="SNV FN",   linewidth=0.5, edgecolor="white")
    ax.bar(x_positions + offsets["indel"], indel_fn_mean, width=bar_width,
           color=indel_colour, label="Indel FN", linewidth=0.5, edgecolor="white")

    ax.set_xticks(x_positions)
    ax.set_xticklabels([f"{v:g}" for v in vafs], rotation=30, ha="right")
    ax.set_xlabel("VAF (%)")
    ax.set_ylabel("Missed truth variants (mean across concentrations)")
    ax.set_ylim(bottom=0)
    ax.legend(frameon=False, loc="upper right")
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_title("Missed variants by type across VAF titration")

    fig.savefig(OUT_DIR / "fn_analysis.pdf")
    fig.savefig(OUT_DIR / "fn_analysis.png", dpi=DPI)
    plt.close(fig)
    print("Saved fn_analysis")


# ── Compare figures ────────────────────────────────────────────────────────────

def _avg_by_vaf(rows, metric):
    """Average a metric across all concentrations at each VAF level.

    Returns a dict of {vaf: mean_value} for VAF > 0.
    """
    from collections import defaultdict
    buckets = defaultdict(list)
    for r in rows:
        if r["vaf"] > 0:
            buckets[r["vaf"]].append(r[metric])
    return {vaf: np.mean(vals) for vaf, vals in buckets.items()}


def _plot_compare_metric(datasets, metric, ylabel, title, out_stem, y_max=None):
    """Line plot comparing a metric across read depths, averaged over concentrations.

    datasets: list of (label, rows) tuples in the order they should appear.
    """
    vafs = sorted({r["vaf"] for _, rows in datasets for r in rows if r["vaf"] > 0})
    x_positions = np.arange(len(vafs))

    fig, ax = plt.subplots(figsize=(COL_W, COL_W * 0.85))

    for (label, rows), colour in zip(datasets, DEPTH_COLOURS):
        avg = _avg_by_vaf(rows, metric)
        yvals = [avg.get(v, 0.0) for v in vafs]
        ax.plot(x_positions, yvals, marker="o", color=colour,
                markersize=4, linewidth=1.0, label=label)

    ax.set_xticks(x_positions)
    ax.set_xticklabels([f"{v:g}" for v in vafs], rotation=30, ha="right")
    ax.set_xlabel("VAF (%)")
    ax.set_ylabel(ylabel)
    if y_max is not None:
        ax.set_ylim(0, y_max)
    else:
        ax.set_ylim(bottom=0)
    ax.yaxis.set_major_formatter(ticker.PercentFormatter(xmax=1, decimals=0))
    ax.legend(title="Read depth", frameon=False, loc="upper left")
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_title(title)

    fig.savefig(OUT_DIR / f"{out_stem}.pdf")
    fig.savefig(OUT_DIR / f"{out_stem}.png", dpi=DPI)
    plt.close(fig)
    print(f"Saved {out_stem}")


def plot_compare(datasets):
    """Generate all read-depth comparison figures.

    datasets: list of (label, rows) tuples, e.g. [("250k", rows), ("500k", rows), ...]
    """
    _plot_compare_metric(
        datasets, "sensitivity",
        ylabel="Sensitivity", title="Sensitivity by read depth (mean across concentrations)",
        out_stem="compare_sensitivity_vs_vaf",
    )
    _plot_compare_metric(
        datasets, "snv_sensitivity",
        ylabel="Sensitivity", title="SNV sensitivity by read depth",
        out_stem="compare_sensitivity_snv_vs_vaf",
    )
    _plot_compare_metric(
        datasets, "indel_sensitivity",
        ylabel="Sensitivity", title="Indel sensitivity by read depth",
        out_stem="compare_sensitivity_indel_vs_vaf",
    )
    _plot_compare_metric(
        datasets, "precision",
        ylabel="Precision", title="Precision by read depth (mean across concentrations)",
        out_stem="compare_precision_vs_vaf",
        y_max=1.0,
    )


# ── Entry point ────────────────────────────────────────────────────────────────

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--results", type=Path, default=_DEFAULT_RESULTS,
                        help="Results TSV for single-run mode "
                             "(default: benchmarking/results/tables_v6a/titration_results_2mreads.tsv)")
    parser.add_argument("--compare", nargs="+", metavar="LABEL:PATH",
                        help="Compare mode: one or more label:path pairs, e.g. "
                             "250k:titration_results_250kreads.tsv. "
                             "Paths relative to benchmarking/results/tables/ if not absolute.")
    return parser.parse_args()


def main():
    args = parse_args()
    tables_dir = REPO / "benchmarking/results/tables"

    if args.compare:
        datasets = []
        for spec in args.compare:
            if ":" not in spec:
                raise SystemExit(f"--compare entries must be LABEL:PATH, got: {spec!r}")
            label, path_str = spec.split(":", 1)
            p = Path(path_str)
            if not p.is_absolute():
                p = tables_dir / p
            rows = load_results(p)
            if not rows:
                print(f"WARNING: no usable rows in {p}; skipping {label!r}")
                continue
            datasets.append((label, rows))
        if not datasets:
            raise SystemExit("No usable datasets for comparison.")
        plot_compare(datasets)
    else:
        rows = load_results(args.results)
        plot_sensitivity(rows)
        plot_precision_recall(rows)
        plot_molecule_stats(rows)
        plot_sensitivity_by_type(rows)
        plot_fp_by_vaf(rows)
        plot_stage_breakdown(rows)
        plot_confusion_matrix(rows)
        plot_fn_analysis(rows)

    print(f"All figures written to {OUT_DIR}")


if __name__ == "__main__":
    main()
