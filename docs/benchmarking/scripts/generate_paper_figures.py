"""
Generate publication-quality benchmark figures for the kam paper.

Reads from:
  docs/benchmarking/summary/sensitivity_by_vaf.csv
  docs/benchmarking/summary/all_results.csv

Writes to:
  docs/paper/figures/
"""

import csv
import os
import sys
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, "..", "..", ".."))
SUMMARY_DIR = os.path.join(REPO_ROOT, "docs", "benchmarking", "summary")
FIGURES_DIR = os.path.join(REPO_ROOT, "docs", "paper", "figures")

SENSITIVITY_CSV = os.path.join(SUMMARY_DIR, "sensitivity_by_vaf.csv")
ALL_RESULTS_CSV = os.path.join(SUMMARY_DIR, "all_results.csv")

os.makedirs(FIGURES_DIR, exist_ok=True)

# ---------------------------------------------------------------------------
# Wong (2011) colourblind-friendly palette (8 colours)
# ---------------------------------------------------------------------------

WONG = {
    "black":          "#000000",
    "orange":         "#E69F00",
    "sky_blue":       "#56B4E9",
    "green":          "#009E73",
    "yellow":         "#F0E442",
    "blue":           "#0072B2",
    "vermillion":     "#D55E00",
    "reddish_purple": "#CC79A7",
}

# ---------------------------------------------------------------------------
# Style constants
# ---------------------------------------------------------------------------

LINE_WIDTH   = 1.5
MARKER_SIZE  = 5
FONT_SIZE    = 10
TICK_SIZE    = 9

plt.rcParams.update({
    "font.size":        FONT_SIZE,
    "axes.titlesize":   FONT_SIZE,
    "axes.labelsize":   FONT_SIZE,
    "xtick.labelsize":  TICK_SIZE,
    "ytick.labelsize":  TICK_SIZE,
    "legend.fontsize":  TICK_SIZE,
    "pdf.fonttype":     42,   # embed TrueType fonts
    "ps.fonttype":      42,
    "figure.dpi":       150,
})

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def save(fig, name):
    """Save figure as both PDF and PNG."""
    pdf_path = os.path.join(FIGURES_DIR, name + ".pdf")
    png_path = os.path.join(FIGURES_DIR, name + ".png")
    fig.savefig(pdf_path, bbox_inches="tight")
    fig.savefig(png_path, bbox_inches="tight", dpi=150)
    print(f"  saved {pdf_path}")
    print(f"  saved {png_path}")


def apply_clean_style(ax):
    """Remove box, add only bottom/left spines, no grid."""
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(direction="out")


def percent_formatter(x, _pos):
    """Format axis tick as percentage string (e.g. 0.001 → '0.1%')."""
    return f"{x * 100:g}%"


def read_sensitivity_csv(path):
    """
    Returns dict: (type, mode) → list of (vaf, mean_sensitivity) sorted by vaf.
    """
    data = defaultdict(list)
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            key = (row["type"], row["mode"])
            data[key].append((float(row["vaf"]), float(row["mean_sensitivity"])))
    for key in data:
        data[key].sort()
    return data


def read_all_results_csv(path):
    """
    Returns list of dicts with keys:
      type, vaf, replicate, mode, tp, fp, fn, sensitivity, precision, f1
    """
    rows = []
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            rows.append({
                "type":        row["type"],
                "vaf":         float(row["vaf"]),
                "replicate":   row["replicate"],
                "mode":        row["mode"],
                "tp":          int(row["tp"]),
                "fp":          int(row["fp"]),
                "fn":          int(row["fn"]),
                "sensitivity": float(row["sensitivity"]),
                "precision":   float(row["precision"]),
                "f1":          float(row["f1"]),
            })
    return rows

# ---------------------------------------------------------------------------
# Figure 1 — Sensitivity vs VAF (main, discovery mode, all types)
# ---------------------------------------------------------------------------

def fig_sensitivity_vs_vaf(sens_data):
    """
    Single-column figure (6 × 4 in).
    One line per variant type in discovery mode.
    X axis: VAF (log scale). Y axis: sensitivity (0–1).
    """
    print("Generating Figure 1: sensitivity vs VAF (all types, discovery)...")

    # Types to plot and their display labels
    type_config = [
        ("snv",    "SNV",               WONG["blue"],      "o"),
        ("indel",  "Indel",             WONG["vermillion"],"s"),
        ("sv",     "SV (DEL/DUP/INV)",  WONG["green"],     "^"),
        ("ins",    "Large insertion",   WONG["orange"],    "D"),
        ("novins", "Novel insertion",   WONG["reddish_purple"], "v"),
    ]

    fig, ax = plt.subplots(figsize=(6, 4))

    for vtype, label, colour, marker in type_config:
        key = (vtype, "discovery")
        if key not in sens_data:
            continue
        points = sens_data[key]
        vafs   = [p[0] for p in points]
        sens   = [p[1] for p in points]
        ax.plot(
            vafs, sens,
            label=label,
            color=colour,
            marker=marker,
            linewidth=LINE_WIDTH,
            markersize=MARKER_SIZE,
            markerfacecolor="white",
            markeredgewidth=LINE_WIDTH,
        )

    ax.set_xscale("log")
    ax.set_xlim(4e-4, 0.12)
    ax.set_ylim(-0.05, 1.05)
    ax.set_xlabel("Variant allele frequency (VAF)")
    ax.set_ylabel("Sensitivity")
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(percent_formatter))
    ax.xaxis.set_minor_formatter(ticker.NullFormatter())
    ax.set_xticks([0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1])

    ax.legend(frameon=False, loc="lower right")
    apply_clean_style(ax)
    fig.tight_layout()
    save(fig, "fig1_sensitivity_vs_vaf")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Figure 2 — SV-specific sensitivity (supplementary)
# ---------------------------------------------------------------------------

def fig_sv_sensitivity(sens_data):
    """
    Single-column figure (6 × 4 in).
    SV subtypes only. Uses aggregate 'sv' if per-subtype rows are absent.
    """
    print("Generating Figure 2 (supplementary): SV sensitivity...")

    sv_config = [
        ("sv",     "SV aggregate (DEL/DUP/INV)", WONG["green"],     "^"),
        ("ins",    "Large insertion",             WONG["orange"],    "D"),
        ("novins", "Novel insertion",             WONG["reddish_purple"], "v"),
    ]

    fig, ax = plt.subplots(figsize=(6, 4))
    plotted = False

    for vtype, label, colour, marker in sv_config:
        key = (vtype, "discovery")
        if key not in sens_data:
            continue
        points = sens_data[key]
        vafs   = [p[0] for p in points]
        sens   = [p[1] for p in points]
        ax.plot(
            vafs, sens,
            label=label,
            color=colour,
            marker=marker,
            linewidth=LINE_WIDTH,
            markersize=MARKER_SIZE,
            markerfacecolor="white",
            markeredgewidth=LINE_WIDTH,
        )
        plotted = True

    if not plotted:
        print("  No SV data found — skipping Figure 2.")
        plt.close(fig)
        return

    ax.set_xscale("log")
    ax.set_xlim(4e-4, 0.12)
    ax.set_ylim(-0.05, 1.05)
    ax.set_xlabel("Variant allele frequency (VAF)")
    ax.set_ylabel("Sensitivity")
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(percent_formatter))
    ax.xaxis.set_minor_formatter(ticker.NullFormatter())
    ax.set_xticks([0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1])

    ax.legend(frameon=False, loc="lower right")
    apply_clean_style(ax)
    fig.tight_layout()
    save(fig, "fig2_sv_sensitivity_supplementary")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Figure 3 — Discovery vs Tumour-Informed sensitivity (if TI data available)
# ---------------------------------------------------------------------------

def fig_discovery_vs_ti(sens_data):
    """
    Double-column figure (10 × 4 in).
    Overlaid curves for discovery and tumour-informed modes.
    Only types that have both modes are included.
    """
    print("Generating Figure 3: discovery vs tumour-informed sensitivity...")

    # Collect types that have tumour_informed data
    ti_types = sorted({vtype for (vtype, mode) in sens_data if mode == "tumour_informed"})
    if not ti_types:
        print("  No tumour-informed data found — skipping Figure 3.")
        return

    type_colours = {
        "ins":    WONG["orange"],
        "novins": WONG["reddish_purple"],
        "sv":     WONG["green"],
        "snv":    WONG["blue"],
        "indel":  WONG["vermillion"],
    }
    type_labels = {
        "ins":    "Large insertion",
        "novins": "Novel insertion",
        "sv":     "SV aggregate",
        "snv":    "SNV",
        "indel":  "Indel",
    }

    fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharey=True)
    ax_disc, ax_ti = axes

    for vtype in ti_types:
        colour = type_colours.get(vtype, WONG["black"])
        label  = type_labels.get(vtype, vtype)

        disc_key = (vtype, "discovery")
        ti_key   = (vtype, "tumour_informed")

        if disc_key in sens_data:
            pts  = sens_data[disc_key]
            vafs = [p[0] for p in pts]
            sens = [p[1] for p in pts]
            ax_disc.plot(vafs, sens, label=label, color=colour,
                         marker="o", linewidth=LINE_WIDTH, markersize=MARKER_SIZE,
                         markerfacecolor="white", markeredgewidth=LINE_WIDTH)

        if ti_key in sens_data:
            pts  = sens_data[ti_key]
            vafs = [p[0] for p in pts]
            sens = [p[1] for p in pts]
            ax_ti.plot(vafs, sens, label=label, color=colour,
                       marker="o", linewidth=LINE_WIDTH, markersize=MARKER_SIZE,
                       markerfacecolor="white", markeredgewidth=LINE_WIDTH)

    for ax, title in [(ax_disc, "Discovery mode"), (ax_ti, "Tumour-informed mode")]:
        ax.set_xscale("log")
        ax.set_xlim(4e-4, 0.12)
        ax.set_ylim(-0.05, 1.05)
        ax.set_xlabel("Variant allele frequency (VAF)")
        ax.xaxis.set_major_formatter(ticker.FuncFormatter(percent_formatter))
        ax.xaxis.set_minor_formatter(ticker.NullFormatter())
        ax.set_xticks([0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1])
        ax.set_title(title)
        apply_clean_style(ax)

    ax_disc.set_ylabel("Sensitivity")
    ax_disc.legend(frameon=False, loc="lower right")

    fig.tight_layout()
    save(fig, "fig3_discovery_vs_tumour_informed")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Figure 4 — False-positive rate by VAF (discovery mode)
# ---------------------------------------------------------------------------

def fig_fp_by_vaf(all_rows):
    """
    Single-column figure (6 × 4 in).
    Mean FP count per sample vs VAF, grouped by variant type.
    Discovery mode only.
    """
    print("Generating Figure 4: FP rate by VAF...")

    # Collect discovery rows only
    disc_rows = [r for r in all_rows if r["mode"] == "discovery"]

    # Compute mean FP per (type, vaf)
    fp_sum   = defaultdict(float)
    fp_count = defaultdict(int)
    for r in disc_rows:
        key = (r["type"], r["vaf"])
        fp_sum[key]   += r["fp"]
        fp_count[key] += 1

    fp_mean = {k: fp_sum[k] / fp_count[k] for k in fp_sum}

    # Collect unique sorted VAFs and types
    all_vafs  = sorted({v for (_t, v) in fp_mean})
    all_types = sorted({t for (t, _v) in fp_mean})

    type_config = {
        "snv":    ("SNV",               WONG["blue"],          "o"),
        "indel":  ("Indel",             WONG["vermillion"],    "s"),
        "sv":     ("SV (DEL/DUP/INV)", WONG["green"],          "^"),
        "ins":    ("Large insertion",   WONG["orange"],        "D"),
        "novins": ("Novel insertion",   WONG["reddish_purple"],"v"),
    }

    fig, ax = plt.subplots(figsize=(6, 4))

    for vtype in ["snv", "indel", "sv", "ins", "novins"]:
        if vtype not in all_types:
            continue
        label, colour, marker = type_config.get(vtype, (vtype, WONG["black"], "o"))
        vafs = []
        fps  = []
        for vaf in all_vafs:
            key = (vtype, vaf)
            if key in fp_mean:
                vafs.append(vaf)
                fps.append(fp_mean[key])
        ax.plot(
            vafs, fps,
            label=label,
            color=colour,
            marker=marker,
            linewidth=LINE_WIDTH,
            markersize=MARKER_SIZE,
            markerfacecolor="white",
            markeredgewidth=LINE_WIDTH,
        )

    ax.set_xscale("log")
    ax.set_xlim(4e-4, 0.12)
    ax.set_xlabel("Variant allele frequency (VAF)")
    ax.set_ylabel("Mean false positives per sample")
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(percent_formatter))
    ax.xaxis.set_minor_formatter(ticker.NullFormatter())
    ax.set_xticks([0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1])

    ax.legend(frameon=False, loc="upper right")
    apply_clean_style(ax)
    fig.tight_layout()
    save(fig, "fig4_fp_by_vaf")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    if not os.path.isfile(SENSITIVITY_CSV):
        print(f"Error: cannot find {SENSITIVITY_CSV}", file=sys.stderr)
        sys.exit(1)
    if not os.path.isfile(ALL_RESULTS_CSV):
        print(f"Error: cannot find {ALL_RESULTS_CSV}", file=sys.stderr)
        sys.exit(1)

    print(f"Reading {SENSITIVITY_CSV}...")
    sens_data = read_sensitivity_csv(SENSITIVITY_CSV)
    print(f"Reading {ALL_RESULTS_CSV}...")
    all_rows = read_all_results_csv(ALL_RESULTS_CSV)

    fig_sensitivity_vs_vaf(sens_data)
    fig_sv_sensitivity(sens_data)
    fig_discovery_vs_ti(sens_data)
    fig_fp_by_vaf(all_rows)

    print("Done.")


if __name__ == "__main__":
    main()
