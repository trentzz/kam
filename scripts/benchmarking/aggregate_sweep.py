#!/usr/bin/env python3
"""Aggregate SV parameter sweep results into summary tables and heatmaps.

Reads result TSVs from kmer_sweep/ and threshold_sweep/ directories,
scores each against truth VCFs, and produces:

  1. A summary CSV with columns:
     parameter, value, vaf, sv_type, sensitivity, precision, n_tp, n_fp

  2. Heatmap CSVs for visualisation (or matplotlib plots if available).

  3. Recommended optimal parameters per SV type.

Usage:
    aggregate_sweep.py --sweep-dir DIR --truth-dir DIR [--output-dir DIR]

Arguments:
    --sweep-dir     Root directory containing kmer_sweep/ and/or threshold_sweep/.
    --truth-dir     Directory containing sim_{type}_vaf{tag}_{rep}/ with truth VCFs.
    --output-dir    Directory for output CSVs (default: sweep-dir/summary).
    -h, --help      Show this help message.
"""

import argparse
import csv
import os
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

POSITION_TOLERANCE = 10  # bp window for position-based matching

VAF_TAGS = [
    "0005", "0010", "0015", "0020", "0025", "0030", "0035", "0040",
    "0050", "0060", "0075", "0100", "0125", "0150", "0175", "0200",
    "0250", "0300", "0350", "0400", "0500", "0600", "0700", "0800",
    "1000",
]

# Map VAF tag to float.
TAG_TO_VAF = {tag: int(tag) / 10000.0 for tag in VAF_TAGS}


def load_vcf_positions(path: Path, filter_pass: bool = False) -> List[int]:
    """Return 1-based POS values from a VCF file."""
    positions = []
    if not path.exists():
        return positions
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.split("\t")
            if len(fields) < 8:
                continue
            if filter_pass and fields[6].strip() != "PASS":
                continue
            try:
                positions.append(int(fields[1]))
            except ValueError:
                pass
    return positions


def score_positions(
    truth_positions: List[int],
    called_positions: List[int],
    tol: int = POSITION_TOLERANCE,
) -> Dict:
    """Score by position proximity. Returns dict with tp, fp, fn, sensitivity, precision."""
    tp = 0
    fn = 0
    for t_pos in truth_positions:
        matched = any(abs(c_pos - t_pos) <= tol for c_pos in called_positions)
        if matched:
            tp += 1
        else:
            fn += 1

    fp = 0
    seen = set()
    for c_pos in called_positions:
        near_truth = any(abs(c_pos - t_pos) <= tol for t_pos in truth_positions)
        if not near_truth and c_pos not in seen:
            fp += 1
            seen.add(c_pos)

    sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    return {
        "n_tp": tp,
        "n_fp": fp,
        "n_fn": fn,
        "sensitivity": round(sensitivity, 4),
        "precision": round(precision, 4),
    }


def find_truth_vcf(truth_dir: Path, sv_type: str, tag: str, rep: str) -> Optional[Path]:
    """Locate the truth VCF for a given sample."""
    sim_dir = truth_dir / f"sim_{sv_type}_vaf{tag}_{rep}"
    if not sim_dir.exists():
        return None
    vcfs = list(sim_dir.glob("*.truth.vcf"))
    return vcfs[0] if vcfs else None


def aggregate_kmer_sweep(
    sweep_dir: Path, truth_dir: Path
) -> List[Dict]:
    """Aggregate results from kmer_sweep/ subdirectories."""
    kmer_dir = sweep_dir / "kmer_sweep"
    if not kmer_dir.exists():
        return []

    rows = []
    for k_dir in sorted(kmer_dir.iterdir()):
        if not k_dir.is_dir():
            continue
        match = re.match(r"k(\d+)", k_dir.name)
        if not match:
            continue
        k = int(match.group(1))

        for result_dir in sorted(k_dir.iterdir()):
            if not result_dir.is_dir():
                continue
            # Parse directory name: {sv_type}_vaf{tag}_{rep}
            m = re.match(r"(\w+)_vaf(\d{4})_([ab])", result_dir.name)
            if not m:
                continue
            sv_type, tag, rep = m.group(1), m.group(2), m.group(3)

            vcf_path = result_dir / "variants.vcf"
            if not vcf_path.exists():
                continue

            truth_vcf = find_truth_vcf(truth_dir, sv_type, tag, rep)
            if truth_vcf is None:
                continue

            truth_pos = load_vcf_positions(truth_vcf)
            called_pos = load_vcf_positions(vcf_path, filter_pass=True)
            scores = score_positions(truth_pos, called_pos)

            rows.append({
                "parameter": "kmer_size",
                "value": str(k),
                "vaf": TAG_TO_VAF.get(tag, 0.0),
                "sv_type": sv_type,
                "replicate": rep,
                **scores,
            })

    return rows


def aggregate_threshold_sweep(
    sweep_dir: Path, truth_dir: Path
) -> List[Dict]:
    """Aggregate results from threshold_sweep/ subdirectories."""
    thresh_dir = sweep_dir / "threshold_sweep"
    if not thresh_dir.exists():
        return []

    rows = []
    for combo_dir in sorted(thresh_dir.iterdir()):
        if not combo_dir.is_dir():
            continue
        # Parse: conf{NNN}_alt{N}
        m = re.match(r"conf(\d{3})_alt(\d+)", combo_dir.name)
        if not m:
            continue
        conf = int(m.group(1)) / 100.0
        alt = int(m.group(2))

        for result_dir in sorted(combo_dir.iterdir()):
            if not result_dir.is_dir():
                continue
            m2 = re.match(r"(\w+)_vaf(\d{4})_([ab])", result_dir.name)
            if not m2:
                continue
            sv_type, tag, rep = m2.group(1), m2.group(2), m2.group(3)

            vcf_path = result_dir / "variants.vcf"
            if not vcf_path.exists():
                continue

            truth_vcf = find_truth_vcf(truth_dir, sv_type, tag, rep)
            if truth_vcf is None:
                continue

            truth_pos = load_vcf_positions(truth_vcf)
            called_pos = load_vcf_positions(vcf_path, filter_pass=True)
            scores = score_positions(truth_pos, called_pos)

            rows.append({
                "parameter": f"conf={conf:.2f}_alt={alt}",
                "value": f"{conf:.2f}/{alt}",
                "vaf": TAG_TO_VAF.get(tag, 0.0),
                "sv_type": sv_type,
                "replicate": rep,
                **scores,
            })

    return rows


def build_heatmap_csv(
    rows: List[Dict], output_dir: Path, parameter_name: str
) -> None:
    """Build a pivot CSV suitable for heatmap visualisation.

    Rows = parameter values, columns = VAF levels, cells = mean sensitivity.
    One file per SV type.
    """
    # Group by (sv_type, parameter_value, vaf).
    grouped = defaultdict(lambda: defaultdict(list))
    for row in rows:
        key = (row["sv_type"], row["value"])
        grouped[key][row["vaf"]].append(row["sensitivity"])

    # Collect all sv_types.
    sv_types = sorted({row["sv_type"] for row in rows})
    all_vafs = sorted({row["vaf"] for row in rows})

    for sv_type in sv_types:
        csv_path = output_dir / f"heatmap_{parameter_name}_{sv_type}.csv"
        param_values = sorted(
            {row["value"] for row in rows if row["sv_type"] == sv_type}
        )

        with open(csv_path, "w", newline="") as fh:
            writer = csv.writer(fh)
            header = [parameter_name] + [f"vaf_{v:.4f}" for v in all_vafs]
            writer.writerow(header)

            for pv in param_values:
                row_data = [pv]
                for vaf in all_vafs:
                    vals = grouped[(sv_type, pv)].get(vaf, [])
                    mean_val = sum(vals) / len(vals) if vals else ""
                    row_data.append(f"{mean_val:.4f}" if isinstance(mean_val, float) else "")
                writer.writerow(row_data)

        print(f"Written: {csv_path}")


def try_plot_heatmaps(output_dir: Path, parameter_name: str) -> None:
    """Attempt to generate matplotlib heatmap plots. Silently skip if unavailable."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError:
        print("[INFO] matplotlib/numpy not available. Skipping heatmap plots.")
        return

    csv_files = list(output_dir.glob(f"heatmap_{parameter_name}_*.csv"))
    for csv_path in csv_files:
        sv_type = csv_path.stem.replace(f"heatmap_{parameter_name}_", "")

        with open(csv_path) as fh:
            reader = csv.reader(fh)
            header = next(reader)
            vaf_labels = [h.replace("vaf_", "") for h in header[1:]]
            param_labels = []
            data = []
            for row in reader:
                param_labels.append(row[0])
                data.append([float(v) if v else 0.0 for v in row[1:]])

        if not data:
            continue

        arr = np.array(data)
        fig, ax = plt.subplots(figsize=(max(8, len(vaf_labels) * 0.5), max(4, len(param_labels) * 0.5)))
        im = ax.imshow(arr, aspect="auto", cmap="YlOrRd", vmin=0, vmax=1)
        ax.set_xticks(range(len(vaf_labels)))
        ax.set_xticklabels(vaf_labels, rotation=45, ha="right", fontsize=7)
        ax.set_yticks(range(len(param_labels)))
        ax.set_yticklabels(param_labels, fontsize=8)
        ax.set_xlabel("VAF")
        ax.set_ylabel(parameter_name)
        ax.set_title(f"Sensitivity: {parameter_name} vs VAF ({sv_type})")
        fig.colorbar(im, ax=ax, label="Sensitivity")
        fig.tight_layout()

        for ext in ("pdf", "png"):
            fig.savefig(output_dir / f"heatmap_{parameter_name}_{sv_type}.{ext}", bbox_inches="tight")
        plt.close(fig)
        print(f"Written: heatmap_{parameter_name}_{sv_type}.pdf/.png")


def recommend_parameters(rows: List[Dict]) -> None:
    """Print recommended optimal parameters per SV type."""
    # Group by (sv_type, parameter, value) and compute mean sensitivity.
    grouped = defaultdict(lambda: defaultdict(list))
    for row in rows:
        key = (row["sv_type"], row["parameter"], row["value"])
        grouped[key].append(row["sensitivity"])

    # For each sv_type, find the parameter value with highest mean sensitivity.
    sv_types = sorted({row["sv_type"] for row in rows})
    param_types = sorted({row["parameter"] for row in rows})

    print("\n=== RECOMMENDED PARAMETERS ===\n")
    for sv_type in sv_types:
        print(f"  {sv_type}:")
        for param in param_types:
            best_value = None
            best_sens = -1.0
            for (st, p, v), sens_list in grouped.items():
                if st == sv_type and p == param:
                    mean_sens = sum(sens_list) / len(sens_list)
                    if mean_sens > best_sens:
                        best_sens = mean_sens
                        best_value = v
            if best_value is not None:
                print(f"    {param}: {best_value} (mean sensitivity: {best_sens:.4f})")
        print()


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Aggregate SV parameter sweep results."
    )
    parser.add_argument(
        "--sweep-dir", required=True,
        help="Root directory containing kmer_sweep/ and/or threshold_sweep/.",
    )
    parser.add_argument(
        "--truth-dir", required=True,
        help="Directory containing sim_{type}_vaf{tag}_{rep}/ with truth VCFs.",
    )
    parser.add_argument(
        "--output-dir", default=None,
        help="Directory for output CSVs (default: sweep-dir/summary).",
    )
    args = parser.parse_args()

    sweep_dir = Path(args.sweep_dir)
    truth_dir = Path(args.truth_dir)
    output_dir = Path(args.output_dir) if args.output_dir else sweep_dir / "summary"
    output_dir.mkdir(parents=True, exist_ok=True)

    all_rows = []

    # Aggregate k-mer sweep.
    kmer_rows = aggregate_kmer_sweep(sweep_dir, truth_dir)
    if kmer_rows:
        print(f"Aggregated {len(kmer_rows)} k-mer sweep results.")
        all_rows.extend(kmer_rows)
        build_heatmap_csv(kmer_rows, output_dir, "kmer_size")
        try_plot_heatmaps(output_dir, "kmer_size")

    # Aggregate threshold sweep.
    thresh_rows = aggregate_threshold_sweep(sweep_dir, truth_dir)
    if thresh_rows:
        print(f"Aggregated {len(thresh_rows)} threshold sweep results.")
        all_rows.extend(thresh_rows)
        build_heatmap_csv(thresh_rows, output_dir, "threshold")
        try_plot_heatmaps(output_dir, "threshold")

    if not all_rows:
        print("[WARN] No sweep results found. Run the sweep scripts first.")
        sys.exit(0)

    # Write combined summary CSV.
    summary_path = output_dir / "sweep_summary.csv"
    fieldnames = [
        "parameter", "value", "vaf", "sv_type", "replicate",
        "sensitivity", "precision", "n_tp", "n_fp", "n_fn",
    ]
    with open(summary_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(all_rows)
    print(f"\nWritten: {summary_path} ({len(all_rows)} rows)")

    # Print recommendations.
    recommend_parameters(all_rows)


if __name__ == "__main__":
    main()
