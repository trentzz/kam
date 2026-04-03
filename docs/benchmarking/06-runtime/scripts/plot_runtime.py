#!/usr/bin/env python3
"""Plot runtime comparison: kam vs alignment pipeline.

Reads both timing CSVs and generates a stacked bar chart showing per-stage
wall-clock times for kam and the alignment pipeline side by side. A secondary
panel shows per-sample wall-clock time as a scatter plot to illustrate variance.

Inputs:
  docs/benchmarking/06-runtime/kam_timings.csv
  docs/benchmarking/06-runtime/alignment_timings.csv

Output:
  docs/benchmarking/06-runtime/figures/runtime_comparison.pdf

CSV column expectations:
  kam_timings.csv:
    sample, wall_time_s, t_assemble_ms, t_index_ms, t_pathfind_ms,
    t_call_ms, t_output_ms, peak_rss_mb

  alignment_timings.csv:
    sample, wall_total_s, t_bwa_s, t_sort_s, t_index_s, t_gatk_s, peak_rss_mb

Run from the repo root:
    python3 docs/benchmarking/06-runtime/scripts/plot_runtime.py
"""

import csv
import sys
from pathlib import Path
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ── Paths ──────────────────────────────────────────────────────────────────────

REPO = Path(__file__).resolve().parents[4]
KAM_CSV   = REPO / "docs/benchmarking/06-runtime/kam_timings.csv"
ALIGN_CSV = REPO / "docs/benchmarking/06-runtime/alignment_timings.csv"
OUT_DIR   = REPO / "docs/benchmarking/06-runtime/figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ── Style ──────────────────────────────────────────────────────────────────────

MM    = 1 / 25.4
COL_W = 84 * MM
DPI   = 600

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

# kam stage colours (Okabe-Ito palette)
KAM_STAGE_COLOURS = {
    "Assemble":  "#0072B2",
    "Index":     "#56B4E9",
    "Pathfind":  "#009E73",
    "Call":      "#E69F00",
    "Output":    "#F0E442",
}

# Alignment stage colours
ALIGN_STAGE_COLOURS = {
    "BWA-MEM2":  "#D55E00",
    "Sort":      "#CC79A7",
    "Index":     "#999999",
    "GATK":      "#7B2D8B",
}


# ── Data loading ──────────────────────────────────────────────────────────────

def _safe_float(val: str) -> float:
    """Return float or 0.0 if the value is empty or non-numeric."""
    try:
        return float(val)
    except (ValueError, TypeError):
        return 0.0


def load_kam_timings(path: Path) -> list:
    """Load kam timings CSV. Returns a list of dicts with numeric fields."""
    if not path.exists():
        print(f"[WARN] Not found: {path}", file=sys.stderr)
        return []
    rows = []
    with open(path) as fh:
        for row in csv.DictReader(fh):
            rows.append({
                "sample":        row["sample"],
                "wall_time_s":   _safe_float(row.get("wall_time_s")),
                "t_assemble_ms": _safe_float(row.get("t_assemble_ms")),
                "t_index_ms":    _safe_float(row.get("t_index_ms")),
                "t_pathfind_ms": _safe_float(row.get("t_pathfind_ms")),
                "t_call_ms":     _safe_float(row.get("t_call_ms")),
                "t_output_ms":   _safe_float(row.get("t_output_ms")),
                "peak_rss_mb":   _safe_float(row.get("peak_rss_mb")),
            })
    return rows


def load_align_timings(path: Path) -> list:
    """Load alignment pipeline timings CSV. Returns a list of dicts."""
    if not path.exists():
        print(f"[WARN] Not found: {path}", file=sys.stderr)
        return []
    rows = []
    with open(path) as fh:
        for row in csv.DictReader(fh):
            rows.append({
                "sample":       row["sample"],
                "wall_total_s": _safe_float(row.get("wall_total_s")),
                "t_bwa_s":      _safe_float(row.get("t_bwa_s")),
                "t_sort_s":     _safe_float(row.get("t_sort_s")),
                "t_index_s":    _safe_float(row.get("t_index_s")),
                "t_gatk_s":     _safe_float(row.get("t_gatk_s")),
                "peak_rss_mb":  _safe_float(row.get("peak_rss_mb")),
            })
    return rows


# ── Aggregate helpers ─────────────────────────────────────────────────────────

def median(values: list) -> float:
    if not values:
        return 0.0
    s = sorted(values)
    n = len(s)
    mid = n // 2
    return s[mid] if n % 2 else (s[mid - 1] + s[mid]) / 2


def kam_stage_medians_s(rows: list) -> dict:
    """Return median per-stage time in seconds for kam."""
    keys = {
        "Assemble":  "t_assemble_ms",
        "Index":     "t_index_ms",
        "Pathfind":  "t_pathfind_ms",
        "Call":      "t_call_ms",
        "Output":    "t_output_ms",
    }
    result = {}
    for label, key in keys.items():
        vals = [r[key] / 1000.0 for r in rows if r[key] > 0]
        result[label] = median(vals)
    return result


def align_stage_medians_s(rows: list) -> dict:
    """Return median per-stage time in seconds for the alignment pipeline."""
    keys = {
        "BWA-MEM2": "t_bwa_s",
        "Sort":     "t_sort_s",
        "Index":    "t_index_s",
        "GATK":     "t_gatk_s",
    }
    result = {}
    for label, key in keys.items():
        vals = [r[key] for r in rows if r[key] > 0]
        result[label] = median(vals)
    return result


# ── Figure: stacked bar chart + per-sample scatter ────────────────────────────

def plot_runtime_comparison(
    kam_rows: list,
    align_rows: list,
    out_dir: Path,
) -> None:
    """Two-panel figure.

    Left: stacked bar chart comparing median per-stage wall-clock time for
          kam (left bar) and the alignment pipeline (right bar).
    Right: per-sample total wall-clock time as a dot plot, one dot per sample.
           kam dots on the left, alignment dots on the right within each pair.
    """
    has_kam   = bool(kam_rows)
    has_align = bool(align_rows)

    if not has_kam and not has_align:
        print(
            "[ERROR] No timing data found. Run time_kam.sh and time_alignment.sh first.",
            file=sys.stderr,
        )
        sys.exit(1)

    with plt.rc_context(RCPARAMS):
        fig, (ax_bar, ax_dot) = plt.subplots(
            1, 2,
            figsize=(COL_W * 2.2, COL_W * 1.0),
            gridspec_kw={"width_ratios": [1.4, 1]},
        )

        # ── Left panel: stacked median-time bars ──────────────────────────────

        bar_x = []
        bar_labels = []
        bar_idx = 0

        # kam bar
        if has_kam:
            kam_stages = kam_stage_medians_s(kam_rows)
            kam_order  = ["Assemble", "Index", "Pathfind", "Call", "Output"]
            bottom = 0.0
            for stage in kam_order:
                h = kam_stages.get(stage, 0.0)
                ax_bar.bar(
                    bar_idx, h, bottom=bottom, width=0.55,
                    color=KAM_STAGE_COLOURS[stage],
                    label=f"kam: {stage}",
                    linewidth=0.3, edgecolor="white",
                )
                bottom += h
            bar_x.append(bar_idx)
            bar_labels.append("kam")
            bar_idx += 1

        # Alignment bar
        if has_align:
            align_stages = align_stage_medians_s(align_rows)
            align_order  = ["BWA-MEM2", "Sort", "Index", "GATK"]
            bottom = 0.0
            for stage in align_order:
                h = align_stages.get(stage, 0.0)
                ax_bar.bar(
                    bar_idx, h, bottom=bottom, width=0.55,
                    color=ALIGN_STAGE_COLOURS[stage],
                    label=f"Alignment: {stage}",
                    linewidth=0.3, edgecolor="white",
                )
                bottom += h
            bar_x.append(bar_idx)
            bar_labels.append("Alignment\npipeline")
            bar_idx += 1

        ax_bar.set_xticks(bar_x)
        ax_bar.set_xticklabels(bar_labels)
        ax_bar.set_ylabel("Wall-clock time (s, median)")
        ax_bar.set_ylim(bottom=0)
        ax_bar.legend(
            frameon=False,
            loc="upper left",
            fontsize=6,
            ncol=1,
        )

        # ── Right panel: per-sample wall-clock scatter ────────────────────────

        jitter_strength = 0.07

        if has_kam:
            rng = np.random.default_rng(seed=42)
            xvals = rng.uniform(-jitter_strength, jitter_strength, size=len(kam_rows))
            yvals = [r["wall_time_s"] for r in kam_rows]
            ax_dot.scatter(
                xvals, yvals,
                color="#0072B2", s=15, alpha=0.7, linewidths=0.3,
                edgecolors="white", label="kam", zorder=3,
            )
            # Median line
            med = median(yvals)
            ax_dot.hlines(
                med, -0.25, 0.25,
                colors="#0072B2", linewidths=1.2, zorder=4,
            )

        if has_align:
            rng = np.random.default_rng(seed=43)
            xvals = 1.0 + rng.uniform(-jitter_strength, jitter_strength, size=len(align_rows))
            yvals = [r["wall_total_s"] for r in align_rows]
            ax_dot.scatter(
                xvals, yvals,
                color="#D55E00", s=15, alpha=0.7, linewidths=0.3,
                edgecolors="white", label="Alignment pipeline", zorder=3,
            )
            med = median(yvals)
            ax_dot.hlines(
                med, 0.75, 1.25,
                colors="#D55E00", linewidths=1.2, zorder=4,
            )

        dot_x = []
        dot_labels = []
        if has_kam:
            dot_x.append(0)
            dot_labels.append("kam")
        if has_align:
            dot_x.append(1)
            dot_labels.append("Alignment\npipeline")

        ax_dot.set_xticks(dot_x)
        ax_dot.set_xticklabels(dot_labels)
        ax_dot.set_xlim(-0.5, 1.5 if has_align else 0.5)
        ax_dot.set_ylim(bottom=0)
        ax_dot.set_ylabel("Wall-clock time (s)")
        ax_dot.legend(frameon=False, loc="upper right", fontsize=6.5)

        fig.tight_layout()

        for ext in ("pdf", "png"):
            out = out_dir / f"runtime_comparison.{ext}"
            fig.savefig(out, bbox_inches="tight")
        plt.close(fig)

    print(f"Written: {out_dir}/runtime_comparison.pdf/.png")


# ── Entry point ───────────────────────────────────────────────────────────────

def main() -> None:
    kam_rows   = load_kam_timings(KAM_CSV)
    align_rows = load_align_timings(ALIGN_CSV)
    plot_runtime_comparison(kam_rows, align_rows, OUT_DIR)


if __name__ == "__main__":
    main()
