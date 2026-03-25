#!/usr/bin/env python3
"""Unified scoring script for all kam benchmark results.

Reads benchmark results from the following directory trees:

  docs/benchmarking/sv/results/sim_*/discovery/variants.vcf
  docs/benchmarking/sv/results/sim_*/tumour_informed/variants.vcf
  docs/benchmarking/sv_new/results/sim_*/discovery/variants.vcf
  docs/benchmarking/sv_new/results/sim_*/tumour_informed/variants.vcf
  docs/benchmarking/snvindel/results/sim_*/discovery/variants.vcf
  docs/benchmarking/snvindel/results/sim_*/tumour_informed/variants.vcf

Each sim_* directory also contains a *.truth.vcf that defines the expected
variants. Scoring is position-based with a ±10 bp tolerance window. This is
necessary because kam reports partial allele representations for large SVs
rather than exact REF/ALT sequences.

The companion variants.tsv is consulted for variant_type metadata when
reported in the summary.

Suite definitions
-----------------
snv       — 5 SNV targets. TP = PASS call within ±10 bp of a truth SNV.
indel     — 5 indel targets. TP = PASS call within ±10 bp of a truth indel.
sv        — 3 targets (DEL, DUP, INV). TP = PASS call near any truth position.
ins       — 1 insertion target. TP = PASS call near the truth insertion site.
invdel    — 1 inversion-deletion target. TP = PASS call near the truth site.
novins    — 1 novel insertion target. TP = PASS call near the truth site.

Outputs
-------
  docs/benchmarking/summary/all_results.csv
      One row per (type, vaf, replicate, mode) with TP, FP, FN, sensitivity,
      precision, and F1.

  docs/benchmarking/summary/sensitivity_by_vaf.csv
      Mean sensitivity per (type, vaf, mode) averaged over replicates.

Run from the repository root:
    python3 docs/benchmarking/scripts/score_all_results.py
"""

import csv
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

# ── Constants ────────────────────────────────────────────────────────────────

REPO_ROOT = Path("docs/benchmarking")

POSITION_TOLERANCE = 10  # bp window for truth-call matching

# Suites and the benchmark root directories they live in.
SUITE_ROOTS: Dict[str, Path] = {
    "snv":    REPO_ROOT / "snvindel",
    "indel":  REPO_ROOT / "snvindel",
    "sv":     REPO_ROOT / "sv",
    "ins":    REPO_ROOT / "sv",
    "invdel": REPO_ROOT / "sv",
    "novins": REPO_ROOT / "sv_new",
}

OUTPUT_DIR = REPO_ROOT / "summary"

ALL_RESULTS_PATH        = OUTPUT_DIR / "all_results.csv"
SENSITIVITY_BY_VAF_PATH = OUTPUT_DIR / "sensitivity_by_vaf.csv"

ALL_RESULTS_COLS = [
    "type", "vaf", "replicate", "mode",
    "tp", "fp", "fn", "sensitivity", "precision", "f1",
]

SENSITIVITY_COLS = [
    "type", "mode", "vaf", "mean_sensitivity", "n_replicates",
]

# Directory name pattern: sim_<type>_vaf<tag>[_<rep>]
# The replicate suffix is optional (some earlier directories omit it).
DIR_PATTERN = re.compile(
    r"^sim_(?P<type>[a-z]+)_vaf(?P<tag>\d+)(?:_(?P<rep>[a-z]))?$"
)


# ── Helpers ──────────────────────────────────────────────────────────────────

def parse_sim_dir(name: str) -> Optional[Tuple[str, str, str]]:
    """Parse a sim directory name into (type, vaf_tag, replicate).

    Returns None if the name does not match the expected pattern.
    """
    m = DIR_PATTERN.match(name)
    if m is None:
        return None
    sv_type = m.group("type")
    tag     = m.group("tag")
    rep     = m.group("rep") or ""
    return sv_type, tag, rep


def tag_to_float(tag: str) -> float:
    """Convert a 4-digit VAF tag (e.g. '0100') to a float (0.01)."""
    return int(tag) / 10000.0


def load_truth_positions(truth_vcf: Path) -> List[int]:
    """Return 1-based POS values from a truth VCF.

    Skips header lines and malformed records. Returns an empty list if the
    file does not exist.
    """
    positions: List[int] = []
    if not truth_vcf.exists():
        return positions
    with open(truth_vcf) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.split("\t")
            if len(fields) < 2:
                continue
            try:
                positions.append(int(fields[1]))
            except ValueError:
                pass
    return positions


def load_pass_positions_vcf(variants_vcf: Path) -> List[int]:
    """Return 1-based POS values for all PASS calls in a variants.vcf file.

    Reads the standard 8-column VCF format produced by kam-call. Only records
    with FILTER == "PASS" are included. Returns an empty list if the file does
    not exist.
    """
    if not variants_vcf.exists():
        return []

    positions: List[int] = []
    with open(variants_vcf) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.split("\t")
            if len(fields) < 7:
                continue
            filt = fields[6].strip()
            if filt != "PASS":
                continue
            try:
                positions.append(int(fields[1]))
            except ValueError:
                pass
    return positions


def score(
    truth_positions: List[int],
    called_positions: List[int],
    tol: int = POSITION_TOLERANCE,
) -> Dict:
    """Compute TP, FP, FN, sensitivity, precision, and F1.

    A truth position is a TP if at least one called position lies within
    ``tol`` bp of it. A called position is a FP if it lies more than ``tol``
    bp from every truth position. FN = number of unmatched truth positions.

    Multiple PASS calls near the same truth position count as one TP.
    """
    tp = 0
    fn = 0
    for t_pos in truth_positions:
        matched = any(abs(c - t_pos) <= tol for c in called_positions)
        if matched:
            tp += 1
        else:
            fn += 1

    fp = 0
    for c_pos in called_positions:
        near_truth = any(abs(c_pos - t) <= tol for t in truth_positions)
        if not near_truth:
            fp += 1

    sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    precision   = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    f1 = (
        2 * precision * sensitivity / (precision + sensitivity)
        if (precision + sensitivity) > 0
        else 0.0
    )
    return {
        "tp":          tp,
        "fp":          fp,
        "fn":          fn,
        "sensitivity": round(sensitivity, 6),
        "precision":   round(precision,   6),
        "f1":          round(f1,           6),
    }


# ── Discovery ────────────────────────────────────────────────────────────────

def collect_results() -> List[Dict]:
    """Walk all known benchmark roots and collect per-sample scored results.

    Discovers sim_* directories dynamically, so new datasets are picked up
    automatically without code changes. Skips directories where no truth VCF
    or no results VCF is present.
    """
    rows: List[Dict] = []

    # Build the set of unique results directories to search.
    results_dirs: List[Path] = []
    seen: set = set()
    for suite_root in SUITE_ROOTS.values():
        d = suite_root / "results"
        if d.is_dir() and d not in seen:
            seen.add(d)
            results_dirs.append(d)

    for results_dir in sorted(results_dirs):
        for sim_dir in sorted(results_dir.iterdir()):
            if not sim_dir.is_dir():
                continue
            parsed = parse_sim_dir(sim_dir.name)
            if parsed is None:
                continue
            sv_type, vaf_tag, rep = parsed

            if sv_type not in SUITE_ROOTS:
                print(
                    f"[SKIP] Unknown suite type '{sv_type}' in {sim_dir}",
                    file=sys.stderr,
                )
                continue

            # Locate the truth VCF (varforge emits it into the sim directory).
            truth_candidates = list(sim_dir.glob("*.truth.vcf"))
            if not truth_candidates:
                print(f"[SKIP] No truth VCF in {sim_dir}", file=sys.stderr)
                continue
            truth_vcf = truth_candidates[0]
            truth_positions = load_truth_positions(truth_vcf)
            if not truth_positions:
                print(
                    f"[WARN] Empty truth VCF: {truth_vcf}",
                    file=sys.stderr,
                )

            vaf_float = tag_to_float(vaf_tag)

            # Score each mode that has a variants.vcf present.
            for mode in ["discovery", "tumour_informed"]:
                vcf_path = sim_dir / mode / "variants.vcf"
                if not vcf_path.exists():
                    continue

                called_positions = load_pass_positions_vcf(vcf_path)
                m = score(truth_positions, called_positions)
                rows.append({
                    "type":      sv_type,
                    "vaf":       vaf_float,
                    "replicate": rep if rep else "?",
                    "mode":      mode,
                    **m,
                })
                print(
                    f"  {sim_dir.name}/{mode}: "
                    f"tp={m['tp']} fp={m['fp']} fn={m['fn']} "
                    f"sens={m['sensitivity']:.4f}",
                    file=sys.stderr,
                )

    return rows


# ── Summary ──────────────────────────────────────────────────────────────────

def compute_sensitivity_by_vaf(rows: List[Dict]) -> List[Dict]:
    """Average sensitivity over replicates for each (type, vaf, mode) group."""
    groups: Dict[Tuple, List[float]] = defaultdict(list)
    for row in rows:
        key = (row["type"], row["vaf"], row["mode"])
        groups[key].append(row["sensitivity"])

    out: List[Dict] = []
    for (sv_type, vaf, mode), sens_list in sorted(groups.items()):
        mean_sens = sum(sens_list) / len(sens_list)
        out.append({
            "type":             sv_type,
            "mode":             mode,
            "vaf":              vaf,
            "mean_sensitivity": round(mean_sens, 6),
            "n_replicates":     len(sens_list),
        })
    return out


def print_summary_table(rows: List[Dict]) -> None:
    """Print a compact sensitivity table to stdout.

    Groups results by (type, mode) and shows mean sensitivity at a set of
    representative VAF levels.
    """
    key_vafs = [0.0010, 0.0025, 0.0050, 0.0100, 0.0200, 0.0500]

    # Accumulate: groups[type][mode][vaf] → list of sensitivity values.
    groups: Dict[str, Dict[str, Dict[float, List[float]]]] = defaultdict(
        lambda: defaultdict(lambda: defaultdict(list))
    )
    for row in rows:
        if row["vaf"] in key_vafs:
            groups[row["type"]][row["mode"]][row["vaf"]].append(
                row["sensitivity"]
            )

    header_vafs = [f"{int(v * 10000):04d}" for v in key_vafs]
    col_header = "  ".join(f"vaf{h}" for h in header_vafs)
    header = f"{'Type':<10} {'Mode':<18} {col_header}"
    print()
    print(header)
    print("-" * len(header))

    suite_order = ["snv", "indel", "sv", "ins", "invdel", "novins"]
    all_types = suite_order + sorted(
        t for t in groups if t not in suite_order
    )

    for sv_type in all_types:
        if sv_type not in groups:
            continue
        for mode in ["discovery", "tumour_informed"]:
            if mode not in groups[sv_type]:
                continue
            vaf_sens = groups[sv_type][mode]
            cells = []
            for v in key_vafs:
                vals = vaf_sens.get(v, [])
                if vals:
                    cells.append(f"{sum(vals)/len(vals):6.3f}")
                else:
                    cells.append(f"{'  -  ':>6}")
            print(f"{sv_type:<10} {mode:<18} " + "  ".join(cells))
    print()


# ── Entry point ───────────────────────────────────────────────────────────────

def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    print("Scanning benchmark results...", file=sys.stderr)
    rows = collect_results()

    if not rows:
        print(
            "[ERROR] No results found. Run from the repository root.",
            file=sys.stderr,
        )
        sys.exit(1)

    # Write all_results.csv.
    with open(ALL_RESULTS_PATH, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=ALL_RESULTS_COLS)
        writer.writeheader()
        writer.writerows(rows)
    print(f"Written: {ALL_RESULTS_PATH} ({len(rows)} rows)")

    # Write sensitivity_by_vaf.csv.
    sens_rows = compute_sensitivity_by_vaf(rows)
    with open(SENSITIVITY_BY_VAF_PATH, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=SENSITIVITY_COLS)
        writer.writeheader()
        writer.writerows(sens_rows)
    print(f"Written: {SENSITIVITY_BY_VAF_PATH} ({len(sens_rows)} rows)")

    # Print summary to stdout.
    print_summary_table(rows)


if __name__ == "__main__":
    main()
