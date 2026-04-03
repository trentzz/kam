#!/usr/bin/env python3
"""Plot sensitivity curves from the varforge SV benchmark suite.

Reads:
  docs/benchmarking/02-sv-core/summary/sv_per_dataset.tsv

Outputs (to docs/benchmarking/results/graphs/):
  varforge_sv_sensitivity.pdf / .png   — discovery sensitivity for sv, ins, invdel

Run from the repo root.
"""

import csv
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = Path(__file__).resolve().parents[4]
PER_DS = REPO / "docs/benchmarking/02-sv-core/summary/sv_per_dataset.tsv"
OUT_DIR = REPO / "docs/benchmarking/results/graphs"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Okabe-Ito palette
BLUE = "#0072B2"
ORANGE = "#E69F00"
GREEN = "#009E73"

MM_PER_INCH = 25.4
SINGLE_COL_INCH = 84 / MM_PER_INCH

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

SV_LABELS = {
    "sv": "DUP+INV",
    "ins": "Insertion",
    "invdel": "INV+DEL",
}
SV_COLOURS = {"sv": BLUE, "ins": ORANGE, "invdel": GREEN}
SV_MARKERS = {"sv": "o", "ins": "s", "invdel": "^"}


def load_per_dataset(path: Path):
    rows = []
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            row["vaf"] = float(row["vaf"])
            row["sensitivity"] = float(row["sensitivity"])
            rows.append(row)
    return rows


def mean_sensitivity(rows, sv_type, mode):
    acc = defaultdict(list)
    for row in rows:
        if row["type"] == sv_type and row["mode"] == mode:
            acc[row["vaf"]].append(row["sensitivity"])
    vafs = sorted(acc.keys())
    means = [sum(acc[v]) / len(acc[v]) for v in vafs]
    return vafs, means


def main() -> None:
    if not PER_DS.exists():
        print(f"[ERROR] Not found: {PER_DS}", file=sys.stderr)
        sys.exit(1)

    rows = load_per_dataset(PER_DS)

    with plt.rc_context(RCPARAMS):
        fig, ax = plt.subplots(figsize=(SINGLE_COL_INCH, SINGLE_COL_INCH * 0.7))

        for sv_type in ["sv", "ins", "invdel"]:
            vafs, sens = mean_sensitivity(rows, sv_type, "discovery")
            x = [v * 100 for v in vafs]
            y = [s * 100 for s in sens]
            ax.plot(x, y,
                    color=SV_COLOURS[sv_type],
                    marker=SV_MARKERS[sv_type],
                    markersize=3,
                    label=SV_LABELS[sv_type])

        ax.set_xscale("log")
        ax.set_xlabel("VAF (%)")
        ax.set_ylabel("Sensitivity (%)")
        ax.set_ylim(0, 105)
        ax.set_xlim(0.04, 12)

        ax.set_xticks([0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0])
        ax.set_xticklabels(["0.05", "0.1", "0.2", "0.5", "1", "2", "5", "10"])

        ax.legend(frameon=False, loc="upper left")

        fig.tight_layout()
        for ext in ("pdf", "png"):
            out = OUT_DIR / f"varforge_sv_sensitivity.{ext}"
            fig.savefig(out, bbox_inches="tight")
        plt.close(fig)
        print(f"Written: {OUT_DIR}/varforge_sv_sensitivity.pdf/.png")


if __name__ == "__main__":
    main()


