#!/usr/bin/env python3
"""Plot sensitivity curves for the sv_new benchmark suite.

Generates per-type sensitivity curves (InvDel, NovIns, Fusion) and a combined
figure that includes the existing sv suite types (DUP+INV, Insertion, InvDel
from the original suite) alongside the new types. Also produces a discovery vs
tumour-informed comparison figure.

Inputs:
  docs/benchmarking/03-sv-extended/summary/sv_new_per_dataset.tsv
  docs/benchmarking/sv/summary/sv_per_dataset.tsv  (optional, for combined figure)

Outputs (to docs/benchmarking/03-sv-extended/figures/):
  sv_new_sensitivity_per_type.pdf   — one panel per new SV type
  sv_new_combined_sensitivity.pdf   — all types including existing sv suite
  sv_new_discovery_vs_ti.pdf        — discovery vs tumour-informed for all new types

Run from the repo root:
    python3 docs/benchmarking/03-sv-extended/scripts/plot_sv_new_results.py
"""

import csv
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ── Paths ──────────────────────────────────────────────────────────────────────

REPO = Path(__file__).resolve().parents[4]
SV_NEW_PER_DS = REPO / "docs/benchmarking/03-sv-extended/summary/sv_new_per_dataset.tsv"
SV_OLD_PER_DS = REPO / "docs/benchmarking/sv/summary/sv_per_dataset.tsv"
OUT_DIR = REPO / "docs/benchmarking/03-sv-extended/figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ── Style (matches graph-style.md rules) ──────────────────────────────────────

MM = 1 / 25.4
COL_W = 84 * MM   # 84 mm single-column width
DPI = 600

# Okabe-Ito palette (colour-blind friendly)
BLUE    = "#0072B2"
ORANGE  = "#E69F00"
GREEN   = "#009E73"
RED     = "#D55E00"
PURPLE  = "#CC79A7"
CYAN    = "#56B4E9"
YELLOW  = "#F0E442"

RCPARAMS = {
    "font.family":        "DejaVu Sans",
    "font.size":          8,
    "axes.labelsize":     8.5,
    "xtick.labelsize":    7.5,
    "ytick.labelsize":    7.5,
    "legend.fontsize":    7.5,
    "axes.spines.top":    False,
    "axes.spines.right":  False,
    "axes.grid":          True,
    "axes.grid.axis":     "y",
    "grid.color":         "#E8E8E8",
    "grid.linewidth":     0.6,
    "lines.linewidth":    1.0,
    "figure.dpi":         150,
    "savefig.dpi":        DPI,
}

# Per-type display configuration.
SV_NEW_TYPES = ["invdel", "novins", "fusion"]

SV_LABELS = {
    "invdel":  "InvDel",
    "novins":  "NovIns",
    "fusion":  "Fusion",
    # Legacy sv suite types (used in the combined figure).
    "sv":      "DUP+INV",
    "ins":     "Insertion",
}

SV_COLOURS = {
    "invdel":  BLUE,
    "novins":  ORANGE,
    "fusion":  GREEN,
    "sv":      RED,
    "ins":     PURPLE,
}

SV_MARKERS = {
    "invdel":  "o",
    "novins":  "s",
    "fusion":  "^",
    "sv":      "D",
    "ins":     "v",
}


# ── Data loading ──────────────────────────────────────────────────────────────

def load_per_dataset(path: Path) -> list:
    """Load a per-dataset TSV. Returns an empty list if the file does not exist."""
    if not path.exists():
        print(f"[WARN] Not found: {path}", file=sys.stderr)
        return []
    rows = []
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            try:
                row["vaf"] = float(row["vaf"])
                row["sensitivity"] = float(row["sensitivity"])
                rows.append(row)
            except (ValueError, KeyError):
                continue
    return rows


def mean_sensitivity_by_vaf(rows: list, sv_type: str, mode: str) -> tuple:
    """Return (vafs, mean_sensitivities) for a given type and mode."""
    acc: dict = defaultdict(list)
    for row in rows:
        if row["type"] == sv_type and row["mode"] == mode:
            acc[row["vaf"]].append(row["sensitivity"])
    vafs = sorted(acc.keys())
    means = [sum(acc[v]) / len(acc[v]) for v in vafs]
    return vafs, means


# ── Figure 1: per-type sensitivity curves ─────────────────────────────────────

def plot_per_type(rows: list, out_dir: Path) -> None:
    """Three-panel figure: one panel per new SV type, discovery mode only."""
    if not rows:
        print("[WARN] No sv_new data; skipping per-type figure.", file=sys.stderr)
        return

    with plt.rc_context(RCPARAMS):
        fig, axes = plt.subplots(
            1, 3,
            figsize=(COL_W * 3.2, COL_W * 0.85),
            sharey=True,
        )

        for ax, sv_type in zip(axes, SV_NEW_TYPES):
            vafs, sens = mean_sensitivity_by_vaf(rows, sv_type, "discovery")
            if not vafs:
                ax.set_title(SV_LABELS[sv_type])
                ax.set_xlabel("VAF (%)")
                continue

            x = [v * 100 for v in vafs]
            y = [s * 100 for s in sens]
            ax.plot(
                x, y,
                color=SV_COLOURS[sv_type],
                marker=SV_MARKERS[sv_type],
                markersize=3,
                label=SV_LABELS[sv_type],
            )
            ax.set_xscale("log")
            ax.set_xlabel("VAF (%)")
            ax.set_ylim(0, 105)
            ax.set_xlim(0.04, 12)
            ax.set_xticks([0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0])
            ax.set_xticklabels(["0.05", "0.1", "0.2", "0.5", "1", "2", "5", "10"])
            ax.set_title(SV_LABELS[sv_type])
            ax.legend(frameon=False, loc="upper left")

        axes[0].set_ylabel("Sensitivity (%)")
        fig.tight_layout()

        for ext in ("pdf", "png"):
            out = out_dir / f"sv_new_sensitivity_per_type.{ext}"
            fig.savefig(out, bbox_inches="tight")
        plt.close(fig)
    print(f"Written: {out_dir}/sv_new_sensitivity_per_type.pdf/.png")


# ── Figure 2: combined sensitivity (new + existing sv suite) ──────────────────

def plot_combined(sv_new_rows: list, sv_old_rows: list, out_dir: Path) -> None:
    """Combined sensitivity curves for all SV types, discovery mode.

    Plots new types (invdel, novins, fusion) and, when available, the existing
    sv suite types (sv, ins) on a single axis.
    """
    if not sv_new_rows:
        print("[WARN] No sv_new data; skipping combined figure.", file=sys.stderr)
        return

    # Build a unified row list with source-tagged types.
    all_rows = list(sv_new_rows)
    all_rows.extend(sv_old_rows)

    # Determine which types have data.
    types_with_data = {row["type"] for row in all_rows}
    plot_order = [t for t in ["sv", "ins", "invdel", "novins", "fusion"]
                  if t in types_with_data]

    with plt.rc_context(RCPARAMS):
        fig, ax = plt.subplots(figsize=(COL_W * 1.4, COL_W * 0.85))

        for sv_type in plot_order:
            vafs, sens = mean_sensitivity_by_vaf(all_rows, sv_type, "discovery")
            if not vafs:
                continue
            x = [v * 100 for v in vafs]
            y = [s * 100 for s in sens]
            ax.plot(
                x, y,
                color=SV_COLOURS[sv_type],
                marker=SV_MARKERS[sv_type],
                markersize=3,
                label=SV_LABELS.get(sv_type, sv_type),
            )

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
            out = out_dir / f"sv_new_combined_sensitivity.{ext}"
            fig.savefig(out, bbox_inches="tight")
        plt.close(fig)
    print(f"Written: {out_dir}/sv_new_combined_sensitivity.pdf/.png")


# ── Figure 3: discovery vs tumour-informed ────────────────────────────────────

def plot_discovery_vs_ti(rows: list, out_dir: Path) -> None:
    """Three-panel figure: discovery (solid) vs tumour-informed (dashed) per type."""
    if not rows:
        print(
            "[WARN] No sv_new data; skipping discovery vs tumour-informed figure.",
            file=sys.stderr,
        )
        return

    with plt.rc_context(RCPARAMS):
        fig, axes = plt.subplots(
            1, 3,
            figsize=(COL_W * 3.2, COL_W * 0.85),
            sharey=True,
        )

        for ax, sv_type in zip(axes, SV_NEW_TYPES):
            colour = SV_COLOURS[sv_type]
            marker = SV_MARKERS[sv_type]

            for mode, linestyle, label_suffix in [
                ("discovery",       "-",  " (discovery)"),
                ("tumour_informed", "--", " (tumour-informed)"),
            ]:
                vafs, sens = mean_sensitivity_by_vaf(rows, sv_type, mode)
                if not vafs:
                    continue
                x = [v * 100 for v in vafs]
                y = [s * 100 for s in sens]
                ax.plot(
                    x, y,
                    color=colour,
                    marker=marker,
                    markersize=3,
                    linestyle=linestyle,
                    label=SV_LABELS[sv_type] + label_suffix,
                )

            ax.set_xscale("log")
            ax.set_xlabel("VAF (%)")
            ax.set_ylim(0, 105)
            ax.set_xlim(0.04, 12)
            ax.set_xticks([0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0])
            ax.set_xticklabels(["0.05", "0.1", "0.2", "0.5", "1", "2", "5", "10"])
            ax.set_title(SV_LABELS[sv_type])
            ax.legend(frameon=False, loc="upper left", fontsize=6.5)

        axes[0].set_ylabel("Sensitivity (%)")
        fig.tight_layout()

        for ext in ("pdf", "png"):
            out = out_dir / f"sv_new_discovery_vs_ti.{ext}"
            fig.savefig(out, bbox_inches="tight")
        plt.close(fig)
    print(f"Written: {out_dir}/sv_new_discovery_vs_ti.pdf/.png")


# ── Entry point ───────────────────────────────────────────────────────────────

def main() -> None:
    sv_new_rows = load_per_dataset(SV_NEW_PER_DS)
    sv_old_rows = load_per_dataset(SV_OLD_PER_DS)

    if not sv_new_rows and not sv_old_rows:
        print(
            "[ERROR] No data found. Run score_sv_new_suite.py first.",
            file=sys.stderr,
        )
        sys.exit(1)

    plot_per_type(sv_new_rows, OUT_DIR)
    plot_combined(sv_new_rows, sv_old_rows, OUT_DIR)
    plot_discovery_vs_ti(sv_new_rows, OUT_DIR)

    print(f"All figures written to {OUT_DIR}")


if __name__ == "__main__":
    main()
