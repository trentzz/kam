#!/usr/bin/env python3
"""Plot sensitivity curves from the varforge SNV/indel benchmark suite.

Reads:
  docs/benchmarking/01-snvindel/summary/snvindel_per_dataset.tsv

Outputs (to docs/benchmarking/results/graphs/):
  varforge_snv_sensitivity.pdf / .png
  varforge_indel_sensitivity.pdf / .png

Run from the repo root.
"""

import csv
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

REPO = Path(__file__).resolve().parents[4]
PER_DS = REPO / "docs/benchmarking/01-snvindel/summary/snvindel_per_dataset.tsv"
OUT_DIR = REPO / "docs/benchmarking/results/graphs"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Okabe-Ito palette
BLUE = "#0072B2"
ORANGE = "#E69F00"
GREEN = "#009E73"
PINK = "#CC79A7"

VAF_LEVELS = [
    0.0005, 0.0010, 0.0015, 0.0020, 0.0025, 0.0030, 0.0035, 0.0040, 0.0050,
    0.0060, 0.0075, 0.0100, 0.0125, 0.0150, 0.0175, 0.0200, 0.0250, 0.0300,
    0.0350, 0.0400, 0.0500, 0.0600, 0.0700, 0.0800, 0.1000,
]

MM_PER_INCH = 25.4
SINGLE_COL_MM = 84
SINGLE_COL_INCH = SINGLE_COL_MM / MM_PER_INCH

RCPARAMS = {
    "font.family": "DejaVu Sans",
    "font.size": 8,
    "axes.labelsize": 8.5,
    "xtick.labelsize": 7.5,
    "ytick.labelsize": 7.5,
    "legend.fontsize": 7.5,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "axes.grid.axis": "y",
    "grid.color": "#E8E8E8",
    "grid.linewidth": 0.6,
    "lines.linewidth": 1.0,
    "figure.dpi": 150,
    "savefig.dpi": 600,
}


def load_per_dataset(path: Path):
    rows = []
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            row["vaf"] = float(row["vaf"])
            row["sensitivity"] = float(row["sensitivity"])
            rows.append(row)
    return rows


def mean_sensitivity(rows, var_type, mode):
    """Return (vaf_levels, mean_sensitivity) averaged over replicates."""
    acc = defaultdict(list)
    for row in rows:
        if row["type"] == var_type and row["mode"] == mode:
            acc[row["vaf"]].append(row["sensitivity"])
    vafs = sorted(acc.keys())
    means = [sum(acc[v]) / len(acc[v]) for v in vafs]
    return vafs, means


def plot_sensitivity_curve(rows, var_type: str, out_stem: str) -> None:
    vafs_disc, sens_disc = mean_sensitivity(rows, var_type, "discovery")
    vafs_ti, sens_ti = mean_sensitivity(rows, var_type, "tumour_informed")

    with plt.rc_context(RCPARAMS):
        fig, ax = plt.subplots(figsize=(SINGLE_COL_INCH, SINGLE_COL_INCH * 0.7))

        x_disc = [v * 100 for v in vafs_disc]
        x_ti = [v * 100 for v in vafs_ti]

        ax.plot(x_disc, [s * 100 for s in sens_disc],
                color=BLUE, marker="o", markersize=3, label="Discovery")
        ax.plot(x_ti, [s * 100 for s in sens_ti],
                color=ORANGE, marker="s", markersize=3, label="Tumour-informed")

        ax.set_xscale("log")
        ax.set_xlabel("VAF (%)")
        ax.set_ylabel("Sensitivity (%)")
        ax.set_ylim(0, 105)
        ax.set_xlim(0.04, 12)

        # Tick labels at readable VAF points
        ax.set_xticks([0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0])
        ax.set_xticklabels(["0.05", "0.1", "0.2", "0.5", "1", "2", "5", "10"])

        ax.legend(frameon=False, loc="upper left")

        fig.tight_layout()
        for ext in ("pdf", "png"):
            out = OUT_DIR / f"{out_stem}.{ext}"
            fig.savefig(out, bbox_inches="tight")
        plt.close(fig)
        print(f"Written: {OUT_DIR / out_stem}.pdf/.png")


def main() -> None:
    if not PER_DS.exists():
        print(f"[ERROR] Not found: {PER_DS}", file=sys.stderr)
        sys.exit(1)

    rows = load_per_dataset(PER_DS)
    plot_sensitivity_curve(rows, "snv", "varforge_snv_sensitivity")
    plot_sensitivity_curve(rows, "indel", "varforge_indel_sensitivity")


if __name__ == "__main__":
    main()
