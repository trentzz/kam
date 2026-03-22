#!/usr/bin/env python3
"""Generate method/explanation figures for the kam paper.

Produces:
  pipeline_overview.pdf   — horizontal pipeline flow diagram with data types
  read_structure.pdf      — Twist UMI read layout and canonical pairing
  debruijn_concept.pdf    — de Bruijn graph showing reference vs variant path

All figures follow graph-style.md rules:
  - DejaVu Sans, 8pt base font
  - 84mm single-column or 174mm double-column width as appropriate
  - 600 DPI raster + PDF
  - Thin lines (0.6–1.0pt)
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.patheffects as pe
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
import numpy as np
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
OUT_DIR = REPO / "docs/paper/figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)

MM = 1 / 25.4
COL_W  = 84  * MM   # single column (inches)
DCOL_W = 174 * MM   # double column (inches)
DPI = 600

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 8,
    "axes.titlesize": 8,
    "axes.labelsize": 8,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "lines.linewidth": 1.0,
    "axes.linewidth": 0.6,
    "figure.dpi": DPI,
    "savefig.dpi": DPI,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.03,
})


# ── Helpers ────────────────────────────────────────────────────────────────────

def _save(fig, stem):
    fig.savefig(OUT_DIR / f"{stem}.pdf")
    fig.savefig(OUT_DIR / f"{stem}.png", dpi=DPI)
    plt.close(fig)
    print(f"Saved {stem}")


def _rounded_box(ax, cx, cy, w, h, label, sublabel=None,
                 facecolor="#dce9f5", edgecolor="#2166ac", fontsize=8):
    """Draw a rounded rectangle centred at (cx, cy) with an optional sublabel."""
    box = FancyBboxPatch((cx - w / 2, cy - h / 2), w, h,
                         boxstyle="round,pad=0.01",
                         linewidth=0.8, edgecolor=edgecolor, facecolor=facecolor,
                         zorder=3)
    ax.add_patch(box)
    if sublabel:
        ax.text(cx, cy + h * 0.12, label, ha="center", va="center",
                fontsize=fontsize, fontweight="bold", zorder=4)
        ax.text(cx, cy - h * 0.2, sublabel, ha="center", va="center",
                fontsize=fontsize - 1.5, color="#555555", zorder=4)
    else:
        ax.text(cx, cy, label, ha="center", va="center",
                fontsize=fontsize, fontweight="bold", zorder=4)


def _arrow(ax, x0, x1, y, label=None, color="#444444"):
    """Draw a horizontal arrow from x0 to x1 at height y, with optional label above."""
    ax.annotate("", xy=(x1, y), xytext=(x0, y),
                arrowprops=dict(arrowstyle="-|>", color=color,
                                lw=0.8, mutation_scale=8),
                zorder=2)
    if label:
        ax.text((x0 + x1) / 2, y + 0.04, label, ha="center", va="bottom",
                fontsize=6, color=color, style="italic", zorder=4)


# ── Figure 1: Pipeline overview ────────────────────────────────────────────────

def plot_pipeline_overview():
    """Horizontal flow diagram of the full kam pipeline."""
    # Stages: name, sublabel (data type out)
    stages = [
        ("FASTQ",     "R1 / R2"),
        ("Parse",     "ReadPair"),
        ("Assemble",  "Molecule"),
        ("Index",     "KmerIndex"),
        ("Pathfind",  "GraphPath"),
        ("Call",      "ScoredPath"),
        ("VCF / TSV", "VariantCall"),
    ]

    n = len(stages)
    box_w = 0.11
    box_h = 0.30
    gap   = 0.045   # gap between box edge and arrow start/end
    y     = 0.5

    fig, ax = plt.subplots(figsize=(DCOL_W, DCOL_W * 0.22))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    # x centres, equally spaced
    xs = np.linspace(0.06, 0.94, n)

    io_colour  = "#2166ac"
    box_colour = "#dce9f5"
    io_face    = "#cfe2f3"

    for i, (name, sublabel) in enumerate(stages):
        is_io = (i == 0 or i == n - 1)
        fc = io_face if is_io else box_colour
        ec = io_colour
        _rounded_box(ax, xs[i], y, box_w, box_h, name, sublabel,
                     facecolor=fc, edgecolor=ec, fontsize=7)

    # Arrows between boxes, labelled with the data type leaving the source stage.
    arrow_labels = [
        "ReadPair", "Molecule", "KmerIndex", "GraphPath", "ScoredPath", "VariantCall",
    ]
    for i in range(n - 1):
        x0 = xs[i]   + box_w / 2 + gap
        x1 = xs[i+1] - box_w / 2 - gap
        _arrow(ax, x0, x1, y, label=None, color="#333333")

    # Data-type labels beneath each arrow.
    for i, lbl in enumerate(arrow_labels):
        x_mid = (xs[i] + box_w / 2 + xs[i + 1] - box_w / 2) / 2
        ax.text(x_mid, y - box_h / 2 - 0.08, lbl, ha="center", va="top",
                fontsize=5.5, color="#555555", style="italic")

    ax.set_title("kam pipeline: six stages from raw FASTQ to variant calls",
                 fontsize=8, pad=4)

    _save(fig, "pipeline_overview")


# ── Figure 2: Read structure ───────────────────────────────────────────────────

def plot_read_structure():
    """Two-panel figure: (a) Twist 5M2S+T read layout, (b) canonical UMI pairing."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(DCOL_W, DCOL_W * 0.30),
                                   gridspec_kw={"width_ratios": [1.2, 1]})

    # ── Panel (a): read layout ─────────────────────────────────────────────────
    seg_colours = {
        "UMI (5 bp)":   "#d01c8b",
        "Skip (2 bp)":  "#b2b2b2",
        "Template":     "#2166ac",
    }
    # Relative widths (5, 2, ~100 → scale template for visual clarity)
    segs = [("UMI (5 bp)", 5), ("Skip (2 bp)", 2), ("Template", 30)]
    total = sum(w for _, w in segs)

    row_y = {"R1": 0.72, "R2": 0.30}
    bar_h = 0.22

    ax1.set_xlim(0, 1)
    ax1.set_ylim(0, 1)
    ax1.axis("off")
    ax1.set_title("(a)  Read structure  (5M2S+T)", fontsize=7.5, loc="left", pad=3)

    for row_label, cy in row_y.items():
        x = 0.05
        ax1.text(0.02, cy, row_label, ha="left", va="center",
                 fontsize=7, fontweight="bold")
        for seg_name, seg_w in segs:
            frac = seg_w / total * 0.88
            colour = seg_colours[seg_name]
            rect = mpatches.FancyBboxPatch((x, cy - bar_h / 2), frac, bar_h,
                                          boxstyle="square,pad=0",
                                          linewidth=0.5, edgecolor="white",
                                          facecolor=colour, zorder=3)
            ax1.add_patch(rect)
            # Label inside segment
            label_short = seg_name.split(" ")[0]
            ax1.text(x + frac / 2, cy, label_short, ha="center", va="center",
                     fontsize=6, color="white", fontweight="bold", zorder=4)
            x += frac

    # Legend
    handles = [mpatches.Patch(facecolor=c, edgecolor="white", label=l)
               for l, c in seg_colours.items()]
    ax1.legend(handles=handles, loc="lower center", fontsize=6,
               frameon=False, ncol=3, bbox_to_anchor=(0.5, -0.05))

    # ── Panel (b): canonical UMI pairing ──────────────────────────────────────
    ax2.set_xlim(0, 1)
    ax2.set_ylim(0, 1)
    ax2.axis("off")
    ax2.set_title("(b)  Canonical UMI pairing", fontsize=7.5, loc="left", pad=3)

    def umi_box(ax, cx, cy, u1, u2, label, colour):
        _rounded_box(ax, cx, cy, 0.38, 0.18,
                     f"{u1} / {u2}", sublabel=label,
                     facecolor=colour, edgecolor="#444444", fontsize=6.5)

    # Forward read pair
    umi_box(ax2, 0.25, 0.78, "ACGTA", "TGCAT", "fwd read pair", "#d1e5f0")
    # Reverse read pair
    umi_box(ax2, 0.25, 0.42, "TGCAT", "ACGTA", "rev read pair", "#fde0ef")
    # Canonical output
    umi_box(ax2, 0.72, 0.60, "ACGTA", "TGCAT", "canonical key", "#e0f5e0")

    # Arrows to canonical
    _arrow(ax2, 0.455, 0.525, 0.78, color="#666666")
    _arrow(ax2, 0.455, 0.525, 0.42, color="#666666")
    # Converge lines
    ax2.annotate("", xy=(0.72, 0.65), xytext=(0.555, 0.78),
                 arrowprops=dict(arrowstyle="-|>", color="#666666",
                                 lw=0.7, mutation_scale=7))
    ax2.annotate("", xy=(0.72, 0.55), xytext=(0.555, 0.42),
                 arrowprops=dict(arrowstyle="-|>", color="#666666",
                                 lw=0.7, mutation_scale=7))

    ax2.text(0.25, 0.12,
             "min(A/B, B/A) lexicographically",
             ha="center", va="center", fontsize=6, color="#555555", style="italic")

    _save(fig, "read_structure")


# ── Figure 3: de Bruijn graph concept ─────────────────────────────────────────

def plot_debruijn_concept():
    """Small de Bruijn graph illustrating a reference path and a SNV variant path."""
    # Manually specify a small k=4 graph.
    # Nodes are k-mers; edges are (k+1)-mers.
    # Reference path: A → B → C → D → E
    # Variant branch:     B → V → D   (SNV causes a different middle k-mer V)

    fig, ax = plt.subplots(figsize=(COL_W, COL_W * 0.85))
    ax.set_xlim(-0.1, 1.1)
    ax.set_ylim(-0.1, 1.1)
    ax.axis("off")
    ax.set_title("de Bruijn graph: reference and variant paths", fontsize=8, pad=4)

    # Node positions
    nodes = {
        "A":  (0.05, 0.50),
        "B":  (0.28, 0.50),
        "C":  (0.55, 0.72),  # reference middle
        "V":  (0.55, 0.28),  # variant middle
        "D":  (0.75, 0.50),
        "E":  (0.97, 0.50),
    }

    # Node labels (representative k-mer sequences)
    labels = {
        "A":  "ACGT",
        "B":  "CGTA",
        "C":  "GTAT",   # reference allele
        "V":  "GTCT",   # variant allele (A→C SNV)
        "D":  "TATG",
        "E":  "ATGC",
    }

    ref_colour  = "#2166ac"
    var_colour  = "#d01c8b"
    node_colour = "#dce9f5"
    node_border = "#2166ac"
    anchor_colour = "#e0f5e0"
    anchor_border = "#1b7837"

    r = 0.065  # node radius

    def draw_node(name, is_anchor=False):
        x, y = nodes[name]
        fc = anchor_colour if is_anchor else node_colour
        ec = anchor_border if is_anchor else node_border
        circle = plt.Circle((x, y), r, facecolor=fc, edgecolor=ec,
                             linewidth=0.8, zorder=3)
        ax.add_patch(circle)
        ax.text(x, y + 0.005, labels[name], ha="center", va="center",
                fontsize=6.5, fontweight="bold", zorder=4,
                color="#1a1a1a")
        # Node name above
        ax.text(x, y + r + 0.035, name, ha="center", va="bottom",
                fontsize=6, color="#444444", zorder=4)

    def draw_edge(n1, n2, colour, style="-", lw=1.0, offset=0.0):
        x0, y0 = nodes[n1]
        x1, y1 = nodes[n2]
        # Shorten by node radius
        dx, dy = x1 - x0, y1 - y0
        length = np.hypot(dx, dy)
        ux, uy = dx / length, dy / length
        sx, sy = x0 + ux * r, y0 + uy * r
        ex, ey = x1 - ux * r, y1 - uy * r
        # Perpendicular offset for overlapping edges
        if offset:
            px, py = -uy * offset, ux * offset
            sx, sy = sx + px, sy + py
            ex, ey = ex + px, ey + py
        ax.annotate("", xy=(ex, ey), xytext=(sx, sy),
                    arrowprops=dict(arrowstyle="-|>", color=colour,
                                    lw=lw, linestyle=style, mutation_scale=9),
                    zorder=2)

    # Draw edges
    # Reference path: A → B → C → D → E
    draw_edge("A", "B", ref_colour, lw=1.1)
    draw_edge("B", "C", ref_colour, lw=1.1)
    draw_edge("C", "D", ref_colour, lw=1.1)
    draw_edge("D", "E", ref_colour, lw=1.1)
    # Variant path: B → V → D
    draw_edge("B", "V", var_colour, lw=1.1)
    draw_edge("V", "D", var_colour, lw=1.1)

    # Draw nodes (anchors A and E have green border)
    for name in ["B", "C", "V", "D"]:
        draw_node(name, is_anchor=False)
    for name in ["A", "E"]:
        draw_node(name, is_anchor=True)

    # Legend
    ref_patch = mpatches.Patch(facecolor=ref_colour, label="Reference path")
    var_patch = mpatches.Patch(facecolor=var_colour, label="Variant path (SNV)")
    anc_patch = mpatches.Patch(facecolor=anchor_colour, edgecolor=anchor_border,
                               linewidth=0.8, label="Anchor k-mer")
    ax.legend(handles=[ref_patch, var_patch, anc_patch],
              loc="lower right", fontsize=6.5, frameon=False)

    # Annotation: the SNV position
    vx, vy = nodes["V"]
    ax.annotate("A→C\nSNV", xy=(vx - 0.01, vy - r - 0.01),
                xytext=(vx - 0.14, vy - 0.22),
                fontsize=6, color=var_colour,
                arrowprops=dict(arrowstyle="->", color=var_colour,
                                lw=0.7, mutation_scale=7))

    _save(fig, "debruijn_concept")


# ── Entry point ────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    plot_pipeline_overview()
    plot_read_structure()
    plot_debruijn_concept()
    print(f"All method figures written to {OUT_DIR}")
