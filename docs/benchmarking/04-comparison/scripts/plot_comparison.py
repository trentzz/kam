"""
plot_comparison.py

Generate publication-quality comparison figures for alignment-based vs kam
variant calling performance.

Reads concordance_summary.csv (and optionally concordance.csv) produced by
build_concordance.py.

Figures produced
----------------
  sensitivity_comparison.pdf  — grouped bar chart of sensitivity by VAF level,
                                 faceted by input amount.
  precision_comparison.pdf    — bar chart of precision / FP rate at 0% VAF.
  concordance_heatmap.pdf     — per-variant per-sample concordance heatmap.
  detection_overlap.pdf       — detection counts per VAF level (stacked bars).

Usage
-----
  python plot_comparison.py [--summary PATH] [--concordance PATH]
                            [--out-dir PATH] [--dpi INT]
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np


# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------

SCRIPT_DIR = Path(__file__).resolve().parent
COMPARISON_DIR = SCRIPT_DIR.parent

DEFAULT_SUMMARY = COMPARISON_DIR / "concordance_summary.csv"
DEFAULT_CONCORDANCE = COMPARISON_DIR / "concordance.csv"
DEFAULT_OUT_DIR = COMPARISON_DIR / "figures"

# Ordered VAF labels for display (excludes 0pc, which is precision-only).
VAF_ORDER = ["0p1pc", "0p25pc", "0p5pc", "1pc", "2pc"]
VAF_LABELS = {
    "0pc":     "0%",
    "0p001pc": "0.001%",
    "0p01pc":  "0.01%",
    "0p1pc":   "0.1%",
    "0p25pc":  "0.25%",
    "0p5pc":   "0.5%",
    "1pc":     "1%",
    "2pc":     "2%",
}
VAF_NUMERIC = {
    "0pc":     0.0,
    "0p001pc": 0.001,
    "0p01pc":  0.01,
    "0p1pc":   0.1,
    "0p25pc":  0.25,
    "0p5pc":   0.5,
    "1pc":     1.0,
    "2pc":     2.0,
}

# Colour palette — colourblind-friendly (Wong 2011).
COLOUR_KAM       = "#0072B2"   # blue
COLOUR_ALIGNMENT = "#E69F00"   # amber
COLOUR_BOTH      = "#009E73"   # green
COLOUR_ALIGN_ONLY = "#56B4E9"  # sky blue
COLOUR_KAM_ONLY   = "#D55E00"  # vermilion
COLOUR_NEITHER    = "#999999"  # grey

# Figure typography.
FONT_SIZE_AXIS  = 12
FONT_SIZE_TITLE = 13
FONT_SIZE_TICK  = 10
FONT_SIZE_LEGEND = 10

# Figure dimensions (inches).
FIG_WIDTH_SINGLE  = 7.0
FIG_HEIGHT_SINGLE = 4.5
FIG_WIDTH_HEATMAP = 12.0


# ---------------------------------------------------------------------------
# Matplotlib style helpers
# ---------------------------------------------------------------------------

def _apply_base_style() -> None:
    """Set global matplotlib style for publication figures."""
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.size": FONT_SIZE_AXIS,
        "axes.titlesize": FONT_SIZE_TITLE,
        "axes.labelsize": FONT_SIZE_AXIS,
        "xtick.labelsize": FONT_SIZE_TICK,
        "ytick.labelsize": FONT_SIZE_TICK,
        "legend.fontsize": FONT_SIZE_LEGEND,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": False,
        "figure.dpi": 150,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "pdf.fonttype": 42,   # embed fonts in PDF
        "ps.fonttype": 42,
    })


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def load_summary(path: Path) -> list[dict]:
    """Load concordance_summary.csv."""
    rows = []
    with path.open(newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            rows.append(row)
    return rows


def load_concordance(path: Path) -> list[dict]:
    """Load concordance.csv."""
    rows = []
    with path.open(newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            rows.append(row)
    return rows


# ---------------------------------------------------------------------------
# Figure 1: Sensitivity comparison bar chart
# ---------------------------------------------------------------------------

def plot_sensitivity_comparison(
    summary_rows: list[dict],
    out_path: Path,
) -> None:
    """Grouped bar chart: kam vs alignment sensitivity, one panel per input amount.

    X axis: VAF level.
    Y axis: sensitivity (0–1).
    Two bars per group: kam (blue) and alignment (amber).
    """
    _apply_base_style()

    # Aggregate sensitivity by (input_ng, vaf_level). The summary CSV does not
    # include input_ng breakdowns directly — those are in concordance.csv. Here
    # we use what we have: overall per-(vaf_level, variant_type) from summary.
    # We aggregate across variant types to get overall sensitivity per VAF level.
    agg: dict[str, dict[str, list[float]]] = defaultdict(lambda: {
        "alignment": [], "kam": []
    })

    for row in summary_rows:
        vaf = row.get("vaf_level", "")
        if vaf not in VAF_ORDER:
            continue
        try:
            align_sens = float(row["alignment_sensitivity"])
            kam_sens   = float(row["kam_sensitivity"])
        except (ValueError, KeyError):
            continue
        agg[vaf]["alignment"].append(align_sens)
        agg[vaf]["kam"].append(kam_sens)

    if not agg:
        print("[WARN] No sensitivity data available for Figure 1. Skipping.", file=sys.stderr)
        return

    vaf_labels_present = [v for v in VAF_ORDER if v in agg]
    x = np.arange(len(vaf_labels_present))
    bar_width = 0.35

    fig, ax = plt.subplots(figsize=(FIG_WIDTH_SINGLE, FIG_HEIGHT_SINGLE))

    align_means = [float(np.mean(agg[v]["alignment"])) for v in vaf_labels_present]
    kam_means   = [float(np.mean(agg[v]["kam"]))       for v in vaf_labels_present]

    bars_align = ax.bar(
        x - bar_width / 2, align_means, bar_width,
        label="Alignment-based", color=COLOUR_ALIGNMENT, edgecolor="white", linewidth=0.5,
    )
    bars_kam = ax.bar(
        x + bar_width / 2, kam_means, bar_width,
        label="kam (tumour-informed)", color=COLOUR_KAM, edgecolor="white", linewidth=0.5,
    )

    ax.set_xlabel("VAF level")
    ax.set_ylabel("Sensitivity")
    ax.set_title("Sensitivity: alignment-based vs kam")
    ax.set_xticks(x)
    ax.set_xticklabels([VAF_LABELS.get(v, v) for v in vaf_labels_present])
    ax.set_ylim(0, 1.1)
    ax.legend(frameon=False)

    # Value labels on bars.
    for bar in (*bars_align, *bars_kam):
        h = bar.get_height()
        if h > 0.02:
            ax.text(
                bar.get_x() + bar.get_width() / 2, h + 0.02,
                f"{h:.2f}", ha="center", va="bottom", fontsize=8,
            )

    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)
    print(f"Figure saved: {out_path}")


# ---------------------------------------------------------------------------
# Figure 2: Precision comparison
# ---------------------------------------------------------------------------

def plot_precision_comparison(
    concordance_rows: list[dict],
    out_path: Path,
) -> None:
    """Bar chart comparing precision (FP rate at 0% VAF).

    Shows the false-positive rate of each method at 0% spiked VAF, where any
    detection is by definition a false positive.
    """
    _apply_base_style()

    # Filter to 0% VAF rows.
    zero_rows = [r for r in concordance_rows if r.get("vaf_level") == "0pc"]

    if not zero_rows:
        print("[WARN] No 0% VAF rows found. Skipping precision figure.", file=sys.stderr)
        return

    total = len(zero_rows)
    align_fp = sum(1 for r in zero_rows if r.get("alignment_detected", "").strip().upper() in ("TRUE", "1", "YES"))
    kam_fp   = sum(1 for r in zero_rows if r.get("kam_detected", "").strip().upper() in ("TRUE", "1", "YES"))

    align_fpr = align_fp / total if total > 0 else 0.0
    kam_fpr   = kam_fp   / total if total > 0 else 0.0

    methods = ["Alignment-based", "kam"]
    fprs    = [align_fpr, kam_fpr]
    colours = [COLOUR_ALIGNMENT, COLOUR_KAM]

    fig, ax = plt.subplots(figsize=(4.5, FIG_HEIGHT_SINGLE))

    bars = ax.bar(methods, fprs, color=colours, edgecolor="white", linewidth=0.5, width=0.5)

    ax.set_ylabel("False positive rate (at 0% VAF)")
    ax.set_title("Precision: false positive rate at 0% VAF")
    ax.set_ylim(0, max(max(fprs) * 1.4, 0.05))

    for bar, fpr in zip(bars, fprs):
        ax.text(
            bar.get_x() + bar.get_width() / 2, fpr + 0.001,
            f"{fpr:.4f}", ha="center", va="bottom", fontsize=10,
        )

    # Annotate total variants checked.
    ax.text(
        0.98, 0.97,
        f"n = {total} variant×sample pairs",
        transform=ax.transAxes,
        ha="right", va="top",
        fontsize=9, color="grey",
    )

    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)
    print(f"Figure saved: {out_path}")


# ---------------------------------------------------------------------------
# Figure 3: Concordance heatmap
# ---------------------------------------------------------------------------

def plot_concordance_heatmap(
    concordance_rows: list[dict],
    out_path: Path,
) -> None:
    """Per-variant × per-sample concordance heatmap.

    Rows: unique truth variants (chromosome:position:ref:alt).
    Columns: unique samples.
    Colour encoding:
      both_detected  — green
      alignment_only — sky blue
      kam_only       — vermilion
      neither        — grey
    """
    _apply_base_style()

    # Collect unique variants and samples (in sorted order for determinism).
    variant_keys: list[tuple[str, int, str, str]] = sorted(
        {(r["chromosome"], int(r["position"]), r["ref"], r["alt"])
         for r in concordance_rows},
        key=lambda v: (v[0], v[1]),
    )
    sample_ids: list[str] = sorted({r["sample_id"] for r in concordance_rows})

    if not variant_keys or not sample_ids:
        print("[WARN] No data for concordance heatmap. Skipping.", file=sys.stderr)
        return

    # Encode concordance categories as integers for the colour map.
    cat_to_int = {
        "both_detected":  0,
        "alignment_only": 1,
        "kam_only":       2,
        "neither":        3,
    }

    # Build the matrix (variants × samples).
    variant_index = {v: i for i, v in enumerate(variant_keys)}
    sample_index  = {s: j for j, s in enumerate(sample_ids)}
    matrix = np.full((len(variant_keys), len(sample_ids)), 3, dtype=int)  # default: neither

    for row in concordance_rows:
        vk = (row["chromosome"], int(row["position"]), row["ref"], row["alt"])
        vi = variant_index.get(vk)
        si = sample_index.get(row["sample_id"])
        if vi is not None and si is not None:
            matrix[vi, si] = cat_to_int.get(row["concordance"], 3)

    # Custom colour map.
    cmap = matplotlib.colors.ListedColormap([
        COLOUR_BOTH,        # 0: both_detected
        COLOUR_ALIGN_ONLY,  # 1: alignment_only
        COLOUR_KAM_ONLY,    # 2: kam_only
        COLOUR_NEITHER,     # 3: neither
    ])
    norm = matplotlib.colors.BoundaryNorm([-0.5, 0.5, 1.5, 2.5, 3.5], cmap.N)

    n_variants = len(variant_keys)
    n_samples  = len(sample_ids)

    # Scale figure height proportional to number of variants.
    fig_height = max(6.0, n_variants * 0.12)
    fig_width  = max(FIG_WIDTH_HEATMAP, n_samples * 0.45)

    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    ax.imshow(matrix, aspect="auto", cmap=cmap, norm=norm, interpolation="none")

    # Sample labels on x axis.
    ax.set_xticks(range(n_samples))
    # Shorten sample IDs for readability.
    short_labels = [_shorten_sample_id(s) for s in sample_ids]
    ax.set_xticklabels(short_labels, rotation=45, ha="right", fontsize=7)

    # Suppress individual variant labels if there are too many.
    if n_variants <= 50:
        ax.set_yticks(range(n_variants))
        ax.set_yticklabels(
            [f"{v[0]}:{v[1]} {v[2]}>{v[3]}" for v in variant_keys],
            fontsize=6,
        )
    else:
        ax.set_yticks([])

    ax.set_xlabel("Sample")
    ax.set_ylabel("Truth variant")
    ax.set_title("Per-variant concordance across samples")

    # Legend.
    legend_patches = [
        mpatches.Patch(color=COLOUR_BOTH,       label="Both detected"),
        mpatches.Patch(color=COLOUR_ALIGN_ONLY,  label="Alignment only"),
        mpatches.Patch(color=COLOUR_KAM_ONLY,    label="kam only"),
        mpatches.Patch(color=COLOUR_NEITHER,     label="Neither"),
    ]
    ax.legend(
        handles=legend_patches,
        loc="upper right",
        bbox_to_anchor=(1.0, 1.0),
        frameon=False,
        fontsize=9,
    )

    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)
    print(f"Figure saved: {out_path}")


def _shorten_sample_id(sample_id: str) -> str:
    """Return a compact label from a sample_id string.

    e.g. TWIST_STDV2_5ng_VAF_1pc_DEDUPED_70bp-targets_results
      → 5ng 1pc
    """
    parts = sample_id.split("_")
    ng_part  = next((p for p in parts if p.endswith("ng")), "")
    vaf_part = ""
    for i, p in enumerate(parts):
        if p == "VAF" and i + 1 < len(parts):
            vaf_part = parts[i + 1]
            break
    if ng_part and vaf_part:
        return f"{ng_part} {VAF_LABELS.get(vaf_part, vaf_part)}"
    return sample_id


# ---------------------------------------------------------------------------
# Figure 4: Detection overlap (stacked bar per VAF level)
# ---------------------------------------------------------------------------

def plot_detection_overlap(
    concordance_rows: list[dict],
    out_path: Path,
) -> None:
    """Stacked bar chart showing detection overlap per VAF level.

    Each bar at a given VAF level shows the proportion of:
      both_detected, alignment_only, kam_only, neither.
    """
    _apply_base_style()

    all_vaf_levels = [
        v for v in (
            ["0p001pc", "0p01pc"] + VAF_ORDER
        )
        if any(r["vaf_level"] == v for r in concordance_rows)
    ]

    if not all_vaf_levels:
        print("[WARN] No data for detection overlap figure. Skipping.", file=sys.stderr)
        return

    # Count concordance categories per VAF level.
    counts: dict[str, dict[str, int]] = {
        v: {"both_detected": 0, "alignment_only": 0, "kam_only": 0, "neither": 0}
        for v in all_vaf_levels
    }
    for row in concordance_rows:
        vaf = row.get("vaf_level", "")
        cat = row.get("concordance", "")
        if vaf in counts and cat in counts[vaf]:
            counts[vaf][cat] += 1

    x = np.arange(len(all_vaf_levels))
    bar_width = 0.6

    fig, ax = plt.subplots(figsize=(FIG_WIDTH_SINGLE, FIG_HEIGHT_SINGLE))

    bottoms = np.zeros(len(all_vaf_levels))
    categories = [
        ("both_detected",  COLOUR_BOTH,       "Both detected"),
        ("alignment_only", COLOUR_ALIGN_ONLY,  "Alignment only"),
        ("kam_only",       COLOUR_KAM_ONLY,    "kam only"),
        ("neither",        COLOUR_NEITHER,     "Neither"),
    ]

    for cat_key, colour, label in categories:
        values = np.array([counts[v][cat_key] for v in all_vaf_levels], dtype=float)
        ax.bar(
            x, values, bar_width,
            bottom=bottoms, color=colour, label=label,
            edgecolor="white", linewidth=0.3,
        )
        bottoms += values

    ax.set_xlabel("VAF level")
    ax.set_ylabel("Number of variant×sample pairs")
    ax.set_title("Detection overlap per VAF level")
    ax.set_xticks(x)
    ax.set_xticklabels([VAF_LABELS.get(v, v) for v in all_vaf_levels])
    ax.legend(frameon=False, loc="upper right")

    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)
    print(f"Figure saved: {out_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Generate comparison figures from concordance data."
    )
    parser.add_argument(
        "--summary",
        type=Path,
        default=DEFAULT_SUMMARY,
        help="Path to concordance_summary.csv (default: %(default)s)",
    )
    parser.add_argument(
        "--concordance",
        type=Path,
        default=DEFAULT_CONCORDANCE,
        help="Path to concordance.csv (default: %(default)s)",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=DEFAULT_OUT_DIR,
        help="Output directory for PDF figures (default: %(default)s)",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=300,
        help="Raster DPI for any raster output (PDF uses vector; default: %(default)s)",
    )
    args = parser.parse_args(argv)

    # ------------------------------------------------------------------
    # Check inputs
    # ------------------------------------------------------------------
    summary_available    = args.summary.exists()
    concordance_available = args.concordance.exists()

    if not summary_available and not concordance_available:
        print(
            "[ERROR] Neither concordance_summary.csv nor concordance.csv found.\n"
            "        Run build_concordance.py first.",
            file=sys.stderr,
        )
        sys.exit(1)

    if not summary_available:
        print(
            f"[WARN] concordance_summary.csv not found at {args.summary}.\n"
            "       Figures that require summary data will be skipped.",
            file=sys.stderr,
        )
    if not concordance_available:
        print(
            f"[WARN] concordance.csv not found at {args.concordance}.\n"
            "       Figures that require per-variant data will be skipped.",
            file=sys.stderr,
        )

    args.out_dir.mkdir(parents=True, exist_ok=True)

    matplotlib.rcParams["savefig.dpi"] = args.dpi

    # ------------------------------------------------------------------
    # Load data
    # ------------------------------------------------------------------
    summary_rows    = load_summary(args.summary)    if summary_available    else []
    concordance_rows = load_concordance(args.concordance) if concordance_available else []

    print(f"Loaded {len(summary_rows)} summary rows, {len(concordance_rows)} concordance rows.")

    # ------------------------------------------------------------------
    # Figure 1: Sensitivity comparison
    # ------------------------------------------------------------------
    if summary_rows:
        plot_sensitivity_comparison(
            summary_rows,
            out_path=args.out_dir / "sensitivity_comparison.pdf",
        )
    else:
        print("[INFO] Skipping sensitivity_comparison.pdf (no summary data).")

    # ------------------------------------------------------------------
    # Figure 2: Precision comparison
    # ------------------------------------------------------------------
    if concordance_rows:
        plot_precision_comparison(
            concordance_rows,
            out_path=args.out_dir / "precision_comparison.pdf",
        )
    else:
        print("[INFO] Skipping precision_comparison.pdf (no concordance data).")

    # ------------------------------------------------------------------
    # Figure 3: Concordance heatmap
    # ------------------------------------------------------------------
    if concordance_rows:
        plot_concordance_heatmap(
            concordance_rows,
            out_path=args.out_dir / "concordance_heatmap.pdf",
        )
    else:
        print("[INFO] Skipping concordance_heatmap.pdf (no concordance data).")

    # ------------------------------------------------------------------
    # Figure 4: Detection overlap
    # ------------------------------------------------------------------
    if concordance_rows:
        plot_detection_overlap(
            concordance_rows,
            out_path=args.out_dir / "detection_overlap.pdf",
        )
    else:
        print("[INFO] Skipping detection_overlap.pdf (no concordance data).")

    print("\nAll figures complete.")


if __name__ == "__main__":
    main()
