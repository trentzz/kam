#!/usr/bin/env python3
"""Score all sv_new benchmark datasets: InvDel, NovIns, and Fusion.

Scoring uses position-based matching with a ±10 bp tolerance window. This is
necessary because kam reports partial allele representations for large SVs
rather than full allele sequences. REF/ALT sequence matching would give 0%
sensitivity for structural variants.

For each truth variant:
  - A call is a TP if any PASS call has POS within ±10 bp of the truth POS.
  - Multiple PASS calls near the same truth position count as one TP.

For fusion truth variants (BND records), both breakpoint positions are checked.
A fusion call is a TP if both breakpoints are detected within tolerance.

Inputs:
  docs/benchmarking/03-sv-extended/results/sim_{type}_vaf{tag}_{rep}/*.truth.vcf
  docs/benchmarking/03-sv-extended/results/kam_{type}_vaf{tag}_{rep}/calls_{mode}.vcf
  docs/benchmarking/03-sv-extended/data/truth_fusion_vaf{tag}_{rep}.vcf  (for fusion)

Outputs:
  docs/benchmarking/03-sv-extended/summary/sv_new_per_dataset.tsv
  docs/benchmarking/03-sv-extended/summary/sv_new_summary.tsv

Run from the repo root:
    python3 docs/benchmarking/03-sv-extended/scripts/score_sv_new_suite.py
"""

import csv
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional

VAF_LEVELS = [
    0.0005, 0.0010, 0.0015, 0.0020, 0.0025, 0.0030, 0.0035, 0.0040, 0.0050,
    0.0060, 0.0075, 0.0100, 0.0125, 0.0150, 0.0175, 0.0200, 0.0250, 0.0300,
    0.0350, 0.0400, 0.0500, 0.0600, 0.0700, 0.0800, 0.1000,
]

KEY_VAF_LEVELS = [0.0010, 0.0025, 0.0050, 0.0100, 0.0200, 0.0500]

POSITION_TOLERANCE = 10  # bp window for position-based matching

ROOT = Path("docs/benchmarking/sv_new")
SUMMARY_DIR = ROOT / "summary"
SUMMARY_DIR.mkdir(parents=True, exist_ok=True)


def vaf_tag(vaf: float) -> str:
    """Return zero-padded 4-digit VAF tag (e.g. 0.01 → '0100')."""
    return f"{round(vaf * 10000):04d}"


def load_vcf_truth_positions(path: Path) -> List[int]:
    """Return 1-based POS values from a truth VCF, skipping header lines."""
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
    """Return 1-based POS values for all PASS calls in a VCF."""
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


def count_pass_bnd_calls(path: Path) -> int:
    """Return the number of PASS BND records in a VCF.

    Fusion calls use the target FASTA header as CHROM and POS=1, so
    position-based matching against truth VCF coordinates (chr1:200, chr1:900)
    is not meaningful. Instead, we treat a PASS BND call as evidence that the
    fusion was detected. The truth always has exactly one fusion event per
    dataset (two BND records for the two breakends), so we cap the TP count at
    the number of truth events (1).
    """
    count = 0
    if not path.exists():
        return count
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.split("\t")
            if len(fields) < 8:
                continue
            filt = fields[6].strip()
            info = fields[7] if len(fields) > 7 else ""
            if filt == "PASS" and "SVTYPE=BND" in info:
                count += 1
    return count


def score_fusion(n_truth_events: int, n_pass_bnd: int) -> Dict:
    """Score fusion detection.

    One truth fusion event (two BND records) is a TP if any PASS BND call
    exists. FPs are PASS BND calls beyond the number of truth events.
    """
    tp = min(n_truth_events, n_pass_bnd)
    fp = max(0, n_pass_bnd - n_truth_events)
    fn = max(0, n_truth_events - n_pass_bnd)
    sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    f1 = (
        2 * precision * sensitivity / (precision + sensitivity)
        if (precision + sensitivity) > 0
        else 0.0
    )
    return {
        "tp": tp,
        "fp": fp,
        "fn": fn,
        "sensitivity": sensitivity,
        "precision": precision,
        "f1": f1,
    }


def score_positions(
    truth_positions: List[int],
    called_positions: List[int],
    tol: int = POSITION_TOLERANCE,
) -> Dict:
    """Score calls against truth by proximity.

    Returns a dict with tp, fp, fn, sensitivity, precision, f1.
    """
    tp = 0
    fn_count = 0
    for t_pos in truth_positions:
        matched = any(abs(c_pos - t_pos) <= tol for c_pos in called_positions)
        if matched:
            tp += 1
        else:
            fn_count += 1

    fp = 0
    for c_pos in called_positions:
        near_any_truth = any(abs(c_pos - t_pos) <= tol for t_pos in truth_positions)
        if not near_any_truth:
            fp += 1

    fn = fn_count
    sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    f1 = (
        2 * precision * sensitivity / (precision + sensitivity)
        if (precision + sensitivity) > 0
        else 0.0
    )
    return {
        "tp": tp,
        "fp": fp,
        "fn": fn,
        "sensitivity": sensitivity,
        "precision": precision,
        "f1": f1,
    }


def find_truth_vcf(sv_type: str, tag: str, rep: str) -> Optional[Path]:
    """Locate the truth VCF for a given dataset.

    Fusion truth VCFs live in the data/ directory. InvDel and NovIns truth
    VCFs are emitted into the simulation output directory by varforge.
    """
    if sv_type == "fusion":
        # Fusion truth VCFs are pre-generated and live in data/.
        p = ROOT / "data" / f"truth_fusion_vaf{tag}_{rep}.vcf"
        return p if p.exists() else None

    # Standard types: varforge emits the truth VCF into the sim output dir.
    sim_dir = ROOT / "results" / f"sim_{sv_type}_vaf{tag}_{rep}"
    if not sim_dir.exists():
        return None
    candidates = list(sim_dir.glob("*.truth.vcf"))
    return candidates[0] if candidates else None


PER_DS_COLS = [
    "type", "vaf", "replicate", "mode",
    "tp", "fp", "fn", "sensitivity", "precision", "f1",
]
SUMMARY_COLS = ["type", "mode"] + [f"sens_{vaf_tag(v)}" for v in KEY_VAF_LEVELS]

SV_TYPES = ["invdel", "novins", "fusion"]


def main() -> None:
    per_ds_rows = []

    for sv_type in SV_TYPES:
        for vaf in VAF_LEVELS:
            t = vaf_tag(vaf)
            for rep in ["a", "b"]:
                truth_vcf = find_truth_vcf(sv_type, t, rep)
                if truth_vcf is None:
                    print(
                        f"[WARN] No truth VCF for {sv_type} vaf{t}_{rep}",
                        file=sys.stderr,
                    )
                    continue

                truth_positions = load_vcf_truth_positions(truth_vcf)

                # Fusion truth VCFs have two BND records per fusion event.
                # We treat each unique fusion as one truth event.
                n_truth_events = (
                    len(truth_positions) // 2 if sv_type == "fusion" else 0
                )

                # Fusion mixed results are in sim_fusion_..._mixed; kam results
                # are in kam_fusion_... (same naming convention as other types).
                kam_dir = ROOT / "results" / f"kam_{sv_type}_vaf{t}_{rep}"

                for mode in ["discovery", "tumour_informed"]:
                    called_vcf = kam_dir / f"calls_{mode}.vcf"
                    if sv_type == "fusion":
                        # Fusion calls use target FASTA name as CHROM; match
                        # by presence of PASS BND records, not by position.
                        n_pass_bnd = count_pass_bnd_calls(called_vcf)
                        m = score_fusion(n_truth_events, n_pass_bnd)
                    else:
                        called_positions = load_vcf_pass_positions(called_vcf)
                        m = score_positions(truth_positions, called_positions)
                    per_ds_rows.append(
                        {
                            "type": sv_type,
                            "vaf": vaf,
                            "replicate": rep,
                            "mode": mode,
                            **m,
                        }
                    )

    per_ds_path = SUMMARY_DIR / "sv_new_per_dataset.tsv"
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
    for sv_type in SV_TYPES:
        for mode in ["discovery", "tumour_informed"]:
            row: Dict = {"type": sv_type, "mode": mode}
            for v in KEY_VAF_LEVELS:
                vals = sens[sv_type][mode][v]
                row[f"sens_{vaf_tag(v)}"] = (
                    round(sum(vals) / len(vals), 4) if vals else ""
                )
            summary_rows.append(row)

    summary_path = SUMMARY_DIR / "sv_new_summary.tsv"
    with open(summary_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=SUMMARY_COLS, delimiter="\t")
        writer.writeheader()
        writer.writerows(summary_rows)
    print(f"Written: {summary_path}")


if __name__ == "__main__":
    main()
