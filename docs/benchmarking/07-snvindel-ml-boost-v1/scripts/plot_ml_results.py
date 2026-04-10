#!/usr/bin/env python3
"""Generate benchmark figures for the ML boost v1 evaluation.

Compares three conditions:
  - baseline:  no ML, no tumour-informed filter (from 01-snvindel 2M-read run)
  - ti:        tumour-informed filter only
  - ml_ti:     ML model (twist-duplex-v1) + tumour-informed filter

Usage:
  python3 plot_ml_results.py \\
    --baseline ../../01-snvindel/summary/titration_results_2mreads.tsv \\
    --ti        ../results/titration_2mreads_ti.tsv \\
    --ml-ti     ../results/titration_2mreads_ml_twist_duplex_ti.tsv

Produces (all saved to docs/paper/figures/):
  ml_sensitivity_vs_vaf.pdf        — three-condition sensitivity comparison
  ml_snv_sensitivity_vs_vaf.pdf    — SNV-specific sensitivity
  ml_indel_sensitivity_vs_vaf.pdf  — indel-specific sensitivity
  ml_ml_sensitivity_vs_vaf.pdf     — ML-filtered sensitivity vs standard PASS
  ml_filter_effect.pdf             — ML_PASS vs ML_FILTER call counts at each VAF
  ml_ml_prob_distribution.pdf      — ml_prob distribution for TP vs FN calls (ML+TI only)

All figures follow graph-style.md rules:
  - DejaVu Sans, 8pt base font
  - 84mm (single-column) or 168mm (double-column) width
  - 600 DPI raster + PDF
  - Thin lines (1.0–1.2pt)
  - No dual y-axes
  - Direct labels where possible
"""

import argparse
import csv
import os
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

# ── Paths ─────────────────────────────────────────────────────────────────────
REPO = Path(__file__).resolve().parents[4]
_DEFAULT_BASELINE = REPO / "docs/benchmarking/01-snvindel/summary/titration_results_2mreads.tsv"
_DEFAULT_TI       = REPO / "docs/benchmarking/07-snvindel-ml-boost-v1/results/titration_2mreads_ti.tsv"
_DEFAULT_ML_TI    = REPO / "docs/benchmarking/07-snvindel-ml-boost-v1/results/titration_2mreads_ml_twist_duplex_ti.tsv"
PER_SAMPLE_DIR    = REPO / "docs/benchmarking/07-snvindel-ml-boost-v1/results/per_sample"
OUT_DIR           = REPO / "docs/paper/figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ── Style ─────────────────────────────────────────────────────────────────────
MM = 1 / 25.4
COL_W = 84 * MM
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

# Condition colours.
COND_COLOURS = {
    "baseline": "#aaaaaa",
    "ti":       "#2166ac",
    "ml_ti":    "#d01c8b",
}
COND_LABELS = {
    "baseline": "Baseline (no ML, no TI)",
    "ti":       "TI only",
    "ml_ti":    "ML + TI",
}
COND_MARKERS = {
    "baseline": "s",
    "ti":       "o",
    "ml_ti":    "^",
}

# Concentration colours (matching 01-snvindel plot_results.py).
NG_COLOURS = {"5ng": "#2166ac", "15ng": "#4dac26", "30ng": "#d01c8b"}
NG_MARKERS  = {"5ng": "o", "15ng": "s", "30ng": "^"}

# ── Load data ─────────────────────────────────────────────────────────────────
def load_results(path):
    """Load a titration results TSV. Skips failed rows."""
    rows = []
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            if not row.get("molecules"):
                continue
            try:
                r = {
                    "sample": row["sample"],
                    "ng": row["ng"],
                    "vaf": float(row["vaf"]),
                    "molecules": int(row["molecules"]),
                    "duplex": int(row["duplex"]),
                    "duplex_pct": float(row["duplex_pct"]),
                    "variants_called": int(row["variants_called"]),
                    "tp": int(row["tp"]), "fp": int(row["fp"]),
                    "fn": int(row["fn"]), "tn": int(row["tn"]),
                    "sensitivity": float(row["sensitivity"]),
                    "precision": float(row["precision"]),
                    "f1": float(row["f1"]),
                    "snv_tp": int(row["snv_tp"]),
                    "snv_fp": int(row["snv_fp"]),
                    "snv_fn": int(row["snv_fn"]),
                    "snv_sensitivity": float(row["snv_sensitivity"]),
                    "snv_precision": float(row["snv_precision"]),
                    "indel_tp": int(row["indel_tp"]),
                    "indel_fp": int(row["indel_fp"]),
                    "indel_fn": int(row["indel_fn"]),
                    "indel_sensitivity": float(row["indel_sensitivity"]),
                    "indel_precision": float(row["indel_precision"]),
                }
                # ML-specific columns (present only in ml_ti run).
                def _safe_float(val):
                    try:
                        return float(val)
                    except (TypeError, ValueError):
                        return None

                def _safe_int(val):
                    try:
                        return int(val)
                    except (TypeError, ValueError):
                        return None

                r["ml_tp"]  = _safe_int(row.get("ml_tp"))
                r["ml_fp"]  = _safe_int(row.get("ml_fp"))
                r["ml_fn"]  = _safe_int(row.get("ml_fn"))
                r["ml_sensitivity"]       = _safe_float(row.get("ml_sensitivity"))
                r["ml_precision"]         = _safe_float(row.get("ml_precision"))
                r["ml_f1"]                = _safe_float(row.get("ml_f1"))
                r["ml_snv_sensitivity"]   = _safe_float(row.get("ml_snv_sensitivity"))
                r["ml_indel_sensitivity"] = _safe_float(row.get("ml_indel_sensitivity"))
                r["ml_snv_tp"]            = _safe_int(row.get("ml_snv_tp"))
                r["ml_snv_fp"]            = _safe_int(row.get("ml_snv_fp"))
                r["ml_snv_fn"]            = _safe_int(row.get("ml_snv_fn"))
                r["ml_indel_tp"]          = _safe_int(row.get("ml_indel_tp"))
                r["ml_indel_fp"]          = _safe_int(row.get("ml_indel_fp"))
                r["ml_indel_fn"]          = _safe_int(row.get("ml_indel_fn"))
                rows.append(r)
            except (ValueError, KeyError):
                continue
    return rows


def _avg_by_vaf(rows, metric):
    """Average a metric across concentrations at each positive VAF level."""
    buckets = defaultdict(list)
    for r in rows:
        if r["vaf"] > 0:
            v = r[metric]
            if v is not None:
                buckets[r["vaf"]].append(v)
    return {vaf: np.mean(vals) for vaf, vals in buckets.items() if vals}


# ── Condition comparison: sensitivity line plots ──────────────────────────────

def _plot_condition_comparison(conditions, metric, ylabel, title, out_stem, y_max=None):
    """Line plot comparing a metric across three conditions, averaged over concentrations.

    conditions: list of (label_key, rows) tuples.
    """
    vafs = sorted({r["vaf"] for _, rows in conditions for r in rows if r["vaf"] > 0})
    x_positions = np.arange(len(vafs))

    fig, ax = plt.subplots(figsize=(COL_W, COL_W * 0.85))

    for (key, rows) in conditions:
        avg = _avg_by_vaf(rows, metric)
        yvals = [avg.get(v, 0.0) for v in vafs]
        ax.plot(x_positions, yvals,
                marker=COND_MARKERS[key], color=COND_COLOURS[key],
                markersize=4, linewidth=1.0, label=COND_LABELS[key])

    ax.set_xticks(x_positions)
    ax.set_xticklabels([f"{v:g}" for v in vafs], rotation=30, ha="right")
    ax.set_xlabel("VAF (%)")
    ax.set_ylabel(ylabel)
    if y_max is not None:
        ax.set_ylim(0, y_max)
    else:
        ax.set_ylim(bottom=0)
    ax.yaxis.set_major_formatter(ticker.PercentFormatter(xmax=1, decimals=0))
    ax.legend(frameon=False, loc="upper left", fontsize=6.5)
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_title(title)

    fig.savefig(OUT_DIR / f"{out_stem}.pdf")
    fig.savefig(OUT_DIR / f"{out_stem}.png", dpi=DPI)
    plt.close(fig)
    print(f"Saved {out_stem}")


def plot_sensitivity_comparison(conditions):
    """Three-panel sensitivity comparison: overall, SNV, indel."""
    fig, axes = plt.subplots(1, 3, figsize=(COL_W * 2.0, COL_W * 0.75), sharey=False)

    metrics = [
        ("sensitivity",       "Overall sensitivity"),
        ("snv_sensitivity",   "SNV sensitivity"),
        ("indel_sensitivity", "Indel sensitivity"),
    ]

    vafs = sorted({r["vaf"] for _, rows in conditions for r in rows if r["vaf"] > 0})
    x_positions = np.arange(len(vafs))

    for ax, (metric, subtitle) in zip(axes, metrics):
        for (key, rows) in conditions:
            avg = _avg_by_vaf(rows, metric)
            yvals = [avg.get(v, 0.0) for v in vafs]
            ax.plot(x_positions, yvals,
                    marker=COND_MARKERS[key], color=COND_COLOURS[key],
                    markersize=3.5, linewidth=1.0, label=COND_LABELS[key])
        ax.set_xticks(x_positions)
        ax.set_xticklabels([f"{v:g}" for v in vafs], rotation=35, ha="right")
        ax.set_xlabel("VAF (%)")
        ax.set_ylabel("Sensitivity")
        ax.set_ylim(0, 1.0)
        ax.yaxis.set_major_formatter(ticker.PercentFormatter(xmax=1, decimals=0))
        ax.legend(frameon=False, loc="upper left", fontsize=5.5)
        ax.spines[["top", "right"]].set_visible(False)
        ax.set_title(subtitle)

    fig.suptitle("Sensitivity comparison: baseline vs TI vs ML+TI", fontsize=8, y=1.01)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "ml_sensitivity_comparison.pdf")
    fig.savefig(OUT_DIR / "ml_sensitivity_comparison.png", dpi=DPI)
    plt.close(fig)
    print("Saved ml_sensitivity_comparison")


# ── ML filter effect: standard PASS vs ML-filtered ───────────────────────────

def plot_ml_filter_effect(ml_ti_rows):
    """Compare standard sensitivity vs ML-filtered sensitivity for the ML+TI run.

    Two-panel line plot: overall and by variant type.
    Shows sensitivity using PASS filter vs sensitivity additionally requiring ML_PASS.
    """
    # Check ML columns are present.
    if not any(r.get("ml_sensitivity") is not None for r in ml_ti_rows):
        print("Skipping ml_filter_effect: no ml_sensitivity data")
        return

    vafs = sorted({r["vaf"] for r in ml_ti_rows if r["vaf"] > 0})
    x_positions = np.arange(len(vafs))

    fig, axes = plt.subplots(1, 2, figsize=(COL_W * 1.9, COL_W * 0.75))

    for ax, (std_metric, ml_metric, title) in zip(axes, [
        ("sensitivity",       "ml_sensitivity",       "Overall sensitivity"),
        ("snv_sensitivity",   "ml_snv_sensitivity",   "SNV sensitivity"),
    ]):
        std_avg = _avg_by_vaf(ml_ti_rows, std_metric)
        ml_avg  = _avg_by_vaf(ml_ti_rows, ml_metric)

        ax.plot(x_positions, [std_avg.get(v, 0.0) for v in vafs],
                marker="o", color="#2166ac", markersize=4, linewidth=1.0,
                label="PASS filter")
        ax.plot(x_positions, [ml_avg.get(v, 0.0) for v in vafs],
                marker="^", color="#d01c8b", markersize=4, linewidth=1.0,
                label="PASS + ML\\_PASS filter")

        ax.set_xticks(x_positions)
        ax.set_xticklabels([f"{v:g}" for v in vafs], rotation=35, ha="right")
        ax.set_xlabel("VAF (%)")
        ax.set_ylabel("Sensitivity")
        ax.set_ylim(0, 1.0)
        ax.yaxis.set_major_formatter(ticker.PercentFormatter(xmax=1, decimals=0))
        ax.legend(frameon=False, loc="upper left", fontsize=6)
        ax.spines[["top", "right"]].set_visible(False)
        ax.set_title(title)

    fig.suptitle("Effect of ML filter on sensitivity (ML+TI run)", fontsize=8, y=1.01)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "ml_filter_effect.pdf")
    fig.savefig(OUT_DIR / "ml_filter_effect.png", dpi=DPI)
    plt.close(fig)
    print("Saved ml_filter_effect")


# ── ML call counts: ML_PASS vs ML_FILTER per VAF ─────────────────────────────

def plot_ml_call_counts(ml_ti_rows):
    """Stacked bar chart showing ML_PASS vs ML_FILTER call counts per VAF level.

    At each VAF, shows how many PASS calls are ML_PASS (kept by ML model)
    vs ML_FILTER (rejected by ML model). Averaged across concentrations.
    """
    if not any(r.get("ml_tp") is not None for r in ml_ti_rows):
        print("Skipping ml_call_counts: no ML columns")
        return

    vafs = sorted({r["vaf"] for r in ml_ti_rows if r["vaf"] > 0})
    x_positions = np.arange(len(vafs))

    # At each VAF, ml_tp = calls that are TP and ML_PASS; standard tp = all TP PASS calls.
    # Calls ML rejects = tp - ml_tp (for TP calls) and fp - ml_fp (for FP calls).
    # Since we're in TI mode, fp = 0 always. So:
    #   ML_PASS calls = ml_tp
    #   ML_FILTER calls = tp - ml_tp (these are TP calls that ML rejected = sensitivity loss)
    ml_pass_avg  = defaultdict(list)
    ml_filter_avg = defaultdict(list)
    for r in ml_ti_rows:
        if r["vaf"] > 0 and r.get("tp") is not None and r.get("ml_tp") is not None:
            tp = r["tp"]
            ml_tp = r["ml_tp"]
            ml_pass_avg[r["vaf"]].append(ml_tp)
            ml_filter_avg[r["vaf"]].append(max(0, tp - ml_tp))

    ml_pass_vals  = [np.mean(ml_pass_avg[v]) if v in ml_pass_avg else 0 for v in vafs]
    ml_filter_vals = [np.mean(ml_filter_avg[v]) if v in ml_filter_avg else 0 for v in vafs]

    fig, ax = plt.subplots(figsize=(COL_W, COL_W * 0.85))

    ax.bar(x_positions, ml_pass_vals, width=0.55, color="#1b7837",
           label="ML\\_PASS (confirmed TP)", linewidth=0.5, edgecolor="white")
    ax.bar(x_positions, ml_filter_vals, width=0.55, color="#d6604d",
           bottom=ml_pass_vals, label="ML\\_FILTER (rejected TP)",
           linewidth=0.5, edgecolor="white")

    ax.set_xticks(x_positions)
    ax.set_xticklabels([f"{v:g}" for v in vafs], rotation=30, ha="right")
    ax.set_xlabel("VAF (%)")
    ax.set_ylabel("TP calls (mean across concentrations)")
    ax.set_ylim(bottom=0)
    ax.legend(frameon=False, loc="upper left")
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_title("ML filter outcome per TP call (ML+TI run)")

    fig.savefig(OUT_DIR / "ml_call_counts.pdf")
    fig.savefig(OUT_DIR / "ml_call_counts.png", dpi=DPI)
    plt.close(fig)
    print("Saved ml_call_counts")


# ── ML probability distribution from per-target TSVs ─────────────────────────

def plot_ml_prob_distribution(per_sample_dir):
    """Violin plot of ml_prob for detected (TP) vs missed (FN) variants.

    Loads all per-sample per-target TSV files and collects ml_prob values
    for detected variants (detected=True) vs missed truth variants (detected=False).

    Only includes rows where ml_prob is not NA (i.e. the variant was called
    with an ML score, meaning it was PASS in the standard filter).
    """
    if not per_sample_dir.exists():
        print("Skipping ml_prob_distribution: per_sample_dir not found")
        return

    tp_probs = []   # ml_prob for called variants that are in truth
    fn_probs = []   # ml_prob for missed variants (detected=False, no ml_prob available)
    all_tp_probs = []  # all ml_prob values for TP calls (detected=True, ml_prob != NA)

    for tsv in sorted(per_sample_dir.glob("*.targets.tsv")):
        with open(tsv) as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                detected = row["detected"].lower() in ("true", "1", "yes")
                ml_prob_str = row.get("ml_prob", "NA")
                if detected and ml_prob_str not in ("NA", ".", ""):
                    try:
                        all_tp_probs.append(float(ml_prob_str))
                    except ValueError:
                        pass

    if not all_tp_probs:
        print("Skipping ml_prob_distribution: no ml_prob values found")
        return

    fig, ax = plt.subplots(figsize=(COL_W, COL_W * 0.85))

    # Histogram of ml_prob for detected TP calls.
    ax.hist(all_tp_probs, bins=20, range=(0, 1), color="#2166ac",
            alpha=0.8, linewidth=0.5, edgecolor="white")

    ax.axvline(0.5, color="#d01c8b", linewidth=0.8, linestyle="--",
               label="ML\\_PASS threshold (0.5)")
    ax.set_xlabel("ML probability (class = real variant)")
    ax.set_ylabel("TP calls")
    ax.set_xlim(0, 1)
    ax.set_ylim(bottom=0)
    ax.legend(frameon=False, fontsize=6.5)
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_title("ML probability distribution for detected TP calls")

    fig.savefig(OUT_DIR / "ml_prob_distribution.pdf")
    fig.savefig(OUT_DIR / "ml_prob_distribution.png", dpi=DPI)
    plt.close(fig)
    print("Saved ml_prob_distribution")


# ── Sensitivity at 2% VAF: grouped bar, all conditions ───────────────────────

def plot_2pct_comparison(conditions):
    """Grouped bar chart at 2% VAF: sensitivity by concentration and condition.

    Three groups (5 ng, 15 ng, 30 ng), within each group three bars (one per condition).
    """
    target_vaf = 2.0
    ng_groups = ["5ng", "15ng", "30ng"]
    cond_keys = [k for k, _ in conditions]

    bar_width = 0.22
    offsets = {k: (i - 1) * bar_width for i, k in enumerate(cond_keys)}

    fig, axes = plt.subplots(1, 3, figsize=(COL_W * 1.9, COL_W * 0.85), sharey=False)

    metric_pairs = [
        ("sensitivity",       "Overall"),
        ("snv_sensitivity",   "SNV"),
        ("indel_sensitivity", "Indel"),
    ]

    for ax, (metric, subtitle) in zip(axes, metric_pairs):
        x_positions = np.arange(len(ng_groups))
        for key, rows in conditions:
            by_ng = {r["ng"]: r[metric] for r in rows
                     if abs(r["vaf"] - target_vaf) < 1e-6}
            yvals = [by_ng.get(ng, 0.0) for ng in ng_groups]
            xpos = x_positions + offsets[key]
            bars = ax.bar(xpos, yvals, width=bar_width,
                          color=COND_COLOURS[key], label=COND_LABELS[key],
                          linewidth=0.5, edgecolor="white")
            for bar, val in zip(bars, yvals):
                if val > 0.01:
                    ax.text(bar.get_x() + bar.get_width() / 2,
                            bar.get_height() + 0.005,
                            f"{val:.0%}", ha="center", va="bottom",
                            fontsize=5.0, rotation=90)

        ax.set_xticks(x_positions)
        ax.set_xticklabels(ng_groups)
        ax.set_xlabel("Input concentration")
        ax.set_ylabel("Sensitivity")
        ax.set_ylim(0, 1.0)
        ax.yaxis.set_major_formatter(ticker.PercentFormatter(xmax=1, decimals=0))
        ax.legend(frameon=False, loc="upper left", fontsize=5.5)
        ax.spines[["top", "right"]].set_visible(False)
        ax.set_title(subtitle)

    fig.suptitle("Sensitivity at 2% VAF by condition and concentration", fontsize=8, y=1.01)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "ml_2pct_comparison.pdf")
    fig.savefig(OUT_DIR / "ml_2pct_comparison.png", dpi=DPI)
    plt.close(fig)
    print("Saved ml_2pct_comparison")


# ── Entry point ───────────────────────────────────────────────────────────────

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--baseline", type=Path, default=_DEFAULT_BASELINE,
                        help="Baseline results TSV (no ML, no TI)")
    parser.add_argument("--ti", type=Path, default=_DEFAULT_TI,
                        help="TI-only results TSV")
    parser.add_argument("--ml-ti", type=Path, default=_DEFAULT_ML_TI,
                        help="ML + TI results TSV")
    parser.add_argument("--per-sample-dir", type=Path, default=PER_SAMPLE_DIR,
                        help="Directory containing per-sample per-target TSV files")
    return parser.parse_args()


def main():
    args = parse_args()

    # Load all three conditions.
    conditions = []
    for key, path in [("baseline", args.baseline), ("ti", args.ti), ("ml_ti", args.ml_ti)]:
        if not path.exists():
            print(f"WARNING: {key} results not found at {path}; skipping")
            continue
        rows = load_results(path)
        if not rows:
            print(f"WARNING: no usable rows in {path}; skipping {key!r}")
            continue
        conditions.append((key, rows))
        print(f"Loaded {len(rows)} rows for condition {key!r}")

    if not conditions:
        raise SystemExit("No usable result files found.")

    # Three-panel sensitivity comparison.
    if len(conditions) >= 2:
        plot_sensitivity_comparison(conditions)
        plot_2pct_comparison(conditions)

    # Per-condition single plots.
    for key, rows in conditions:
        _plot_condition_comparison(
            [(key, rows)], "sensitivity",
            ylabel="Sensitivity",
            title=f"Sensitivity vs VAF — {COND_LABELS[key]}",
            out_stem=f"ml_{key}_sensitivity_vs_vaf",
        )

    # ML-specific plots (require ml_ti run).
    ml_ti_rows = next((r for k, r in conditions if k == "ml_ti"), None)
    if ml_ti_rows:
        plot_ml_filter_effect(ml_ti_rows)
        plot_ml_call_counts(ml_ti_rows)
        plot_ml_prob_distribution(args.per_sample_dir)

    print(f"\nAll figures written to {OUT_DIR}")


if __name__ == "__main__":
    main()
