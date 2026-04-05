#!/usr/bin/env python3
"""Score all SNV and indel benchmark datasets.

For each of the 50 SNV and 50 indel datasets (25 VAF levels × 2 replicates),
reads the discovery and tumour-informed VCFs and scores them against the
varforge truth VCF.

Matching is POSITION-BASED with a ±10bp tolerance window. This prevents
false positives from REF/ALT collisions (different positions sharing the same
base substitution) from inflating sensitivity scores.

Outputs:
  docs/benchmarking/01-snvindel/summary/snvindel_per_dataset.tsv
  docs/benchmarking/01-snvindel/summary/snvindel_summary.tsv

Run from the repo root.
"""

import csv
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List

VAF_LEVELS = [
    0.0005, 0.0010, 0.0015, 0.0020, 0.0025, 0.0030, 0.0035, 0.0040, 0.0050,
    0.0060, 0.0075, 0.0100, 0.0125, 0.0150, 0.0175, 0.0200, 0.0250, 0.0300,
    0.0350, 0.0400, 0.0500, 0.0600, 0.0700, 0.0800, 0.1000,
]

KEY_VAF_LEVELS = [0.0010, 0.0025, 0.0050, 0.0100, 0.0200, 0.0500]

POSITION_TOLERANCE = 10  # bp window for position-based matching

ROOT = Path("docs/benchmarking/snvindel")
SUMMARY_DIR = ROOT / "summary"
SUMMARY_DIR.mkdir(parents=True, exist_ok=True)


def vaf_tag(vaf: float) -> str:
    return f"{round(vaf * 10000):04d}"


def load_vcf_truth_positions(path: Path) -> List[int]:
    """Return 1-based POS values from a truth VCF."""
    positions = []
    if not path.exists():
        return positions
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.split("\t")
            if len(fields) < 5:
                continue
            try:
                positions.append(int(fields[1]))
            except ValueError:
                pass
    return positions


def load_vcf_pass_positions(path: Path) -> List[int]:
    """Return 1-based POS values for all PASS calls."""
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
            filt = fields[6].strip()
            if filt == "PASS":
                try:
                    positions.append(int(fields[1]))
                except ValueError:
                    pass
    return positions


def score_positions(truth_positions: List[int], called_positions: List[int],
                    tol: int = POSITION_TOLERANCE) -> Dict:
    """Score by position proximity."""
    tp = 0
    fn_list = []
    for t_pos in truth_positions:
        matched = any(abs(c_pos - t_pos) <= tol for c_pos in called_positions)
        if matched:
            tp += 1
        else:
            fn_list.append(t_pos)

    fp = 0
    matched_called = set()
    for c_pos in called_positions:
        near_any_truth = any(abs(c_pos - t_pos) <= tol for t_pos in truth_positions)
        if not near_any_truth and c_pos not in matched_called:
            fp += 1
            matched_called.add(c_pos)

    fn = len(fn_list)
    sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    f1 = 2 * precision * sensitivity / (precision + sensitivity) if (precision + sensitivity) > 0 else 0.0
    return {"tp": tp, "fp": fp, "fn": fn, "sensitivity": sensitivity,
            "precision": precision, "f1": f1}


PER_DS_COLS = ["type", "vaf", "replicate", "mode", "tp", "fp", "fn",
               "sensitivity", "precision", "f1"]
SUMMARY_COLS = ["type", "mode"] + [f"sens_{vaf_tag(v)}" for v in KEY_VAF_LEVELS]


def main() -> None:
    per_ds_rows = []

    for var_type in ["snv", "indel"]:
        for vaf in VAF_LEVELS:
            t = vaf_tag(vaf)
            for rep in ["a", "b"]:
                sim_dir = ROOT / "results" / f"sim_{var_type}_vaf{t}_{rep}"
                kam_dir = ROOT / "results" / f"kam_{var_type}_vaf{t}_{rep}"

                truth_vcf = next(sim_dir.glob("*.truth.vcf"), None) if sim_dir.exists() else None
                if truth_vcf is None:
                    print(f"[WARN] No truth VCF in {sim_dir}", file=sys.stderr)
                    continue

                truth_positions = load_vcf_truth_positions(truth_vcf)

                for mode in ["discovery", "tumour_informed"]:
                    called_vcf = kam_dir / f"calls_{mode}.vcf"
                    called_positions = load_vcf_pass_positions(called_vcf)
                    m = score_positions(truth_positions, called_positions)
                    per_ds_rows.append({
                        "type": var_type,
                        "vaf": vaf,
                        "replicate": rep,
                        "mode": mode,
                        **m,
                    })

    per_ds_path = SUMMARY_DIR / "snvindel_per_dataset.tsv"
    with open(per_ds_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=PER_DS_COLS, delimiter="\t")
        writer.writeheader()
        writer.writerows(per_ds_rows)
    print(f"Written: {per_ds_path} ({len(per_ds_rows)} rows)")

    # Build summary table: mean sensitivity at key VAF levels.
    sens: Dict = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
    for row in per_ds_rows:
        if row["vaf"] in KEY_VAF_LEVELS:
            sens[row["type"]][row["mode"]][row["vaf"]].append(row["sensitivity"])

    summary_rows = []
    for var_type in ["snv", "indel"]:
        for mode in ["discovery", "tumour_informed"]:
            row = {"type": var_type, "mode": mode}
            for v in KEY_VAF_LEVELS:
                vals = sens[var_type][mode][v]
                row[f"sens_{vaf_tag(v)}"] = round(sum(vals) / len(vals), 4) if vals else ""
            summary_rows.append(row)

    summary_path = SUMMARY_DIR / "snvindel_summary.tsv"
    with open(summary_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=SUMMARY_COLS, delimiter="\t")
        writer.writeheader()
        writer.writerows(summary_rows)
    print(f"Written: {summary_path}")


if __name__ == "__main__":
    main()
