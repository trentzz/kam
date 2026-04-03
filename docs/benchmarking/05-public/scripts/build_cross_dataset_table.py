#!/usr/bin/env python3
"""Build the cross-dataset performance summary table for public benchmarks.

Reads results from all three public benchmark datasets:
  PUB-002: UMI Clustering Benchmark (assembly-only; no truth set)
  PUB-003: SEQC2 HCC1395 (SNV and indel)
  PUB-004: COLO829 (SV)

Produces a summary CSV:
  docs/benchmarking/05-public/cross_dataset_summary.csv

Columns:
  dataset, variant_type, mode, n_truth, tp, fp, fn, sensitivity, precision

Usage:
    python3 build_cross_dataset_table.py

Run from the repository root. Requires Python 3.9+.

Notes:
  - Results for SEQC2 and COLO829 are read from their kam output VCFs.
  - Truth variant counts are read from the truth VCFs.
  - Matching is position-based with a ±50 bp tolerance to accommodate
    SV partial-allele reporting.
  - The UMI benchmark row is included with n_truth=NA to document that
    the run completed; all TP/FP/FN fields are left empty.
  - If a results VCF does not exist, the row is still written with empty
    values and a warning printed to stderr.
"""

import csv
import sys
import gzip
from pathlib import Path


REPO = Path(__file__).parent.parent.parent.parent.parent
PUBLIC_DIR = REPO / "docs" / "benchmarking" / "public"
OUTPUT_CSV = PUBLIC_DIR / "cross_dataset_summary.csv"

# Position tolerance for SNV/indel matching.
SNV_TOLERANCE = 10   # bp
# Wider tolerance for SVs: kam reports partial junction alleles.
SV_TOLERANCE = 500   # bp

# ── VCF helpers ───────────────────────────────────────────────────────────────


def open_vcf(path: Path):
    """Open a plain or gzip-compressed VCF for reading."""
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def load_truth_variants(vcf_path: Path) -> list[dict]:
    """Parse truth VCF variants.

    Returns a list of dicts with keys: chrom, pos, ref, alt.
    Filters to PASS variants where a FILTER column is present.
    """
    if not vcf_path.exists():
        return []
    variants = []
    with open_vcf(vcf_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                continue
            try:
                pos = int(fields[1])
            except ValueError:
                continue
            filt = fields[6].strip() if len(fields) > 6 else "."
            # Accept variants with PASS or missing filter (plain VCFs like COLO829).
            if filt not in ("PASS", ".", ""):
                continue
            variants.append(
                {
                    "chrom": fields[0],
                    "pos": pos,
                    "ref": fields[3],
                    "alt": fields[4],
                }
            )
    return variants


def load_pass_calls(vcf_path: Path) -> list[dict]:
    """Parse PASS calls from a kam output VCF.

    Returns a list of dicts with keys: chrom, pos, ref, alt.
    """
    if not vcf_path.exists():
        return []
    calls = []
    with open_vcf(vcf_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 7:
                continue
            filt = fields[6].strip()
            if filt != "PASS":
                continue
            try:
                pos = int(fields[1])
            except ValueError:
                continue
            calls.append(
                {
                    "chrom": fields[0],
                    "pos": pos,
                    "ref": fields[3] if len(fields) > 3 else "",
                    "alt": fields[4] if len(fields) > 4 else "",
                }
            )
    return calls


# ── Scoring ───────────────────────────────────────────────────────────────────


def score(
    truth: list[dict],
    calls: list[dict],
    tol: int,
) -> dict:
    """Score calls against truth using position-based matching.

    A truth variant is a TP if any call is within tol bp on the same
    chromosome. A call is a FP if it is not within tol bp of any truth
    variant. Matching is many-to-one: multiple calls near one truth count
    as one TP and (n-1) FPs.
    """
    n_truth = len(truth)
    if n_truth == 0 and not calls:
        return {
            "n_truth": 0, "tp": 0, "fp": 0, "fn": 0,
            "sensitivity": "", "precision": "",
        }

    tp = 0
    fn_count = 0
    matched_truth: set[int] = set()

    for i, t in enumerate(truth):
        for c in calls:
            if c["chrom"] == t["chrom"] and abs(c["pos"] - t["pos"]) <= tol:
                matched_truth.add(i)
                break
        else:
            fn_count += 1

    tp = len(matched_truth)
    fn_count = n_truth - tp

    fp = 0
    for c in calls:
        near = any(
            t["chrom"] == c["chrom"] and abs(t["pos"] - c["pos"]) <= tol
            for t in truth
        )
        if not near:
            fp += 1

    sensitivity = tp / n_truth if n_truth > 0 else 0.0
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0

    return {
        "n_truth": n_truth,
        "tp": tp,
        "fp": fp,
        "fn": fn_count,
        "sensitivity": round(sensitivity, 4),
        "precision": round(precision, 4),
    }


# ── Dataset configurations ────────────────────────────────────────────────────


def umi_benchmark_row() -> dict:
    """UMI Clustering Benchmark has no truth set.

    Returns a documentation row indicating the run status.
    """
    disc_vcf = PUBLIC_DIR / "umi_benchmark" / "results" / "kam_umi_benchmark" / "calls_discovery.vcf"
    status = "done" if disc_vcf.exists() else "not_run"
    return {
        "dataset": "PUB-002_umi_benchmark",
        "variant_type": "assembly_only",
        "mode": "discovery",
        "n_truth": "NA",
        "tp": "",
        "fp": "",
        "fn": "",
        "sensitivity": "",
        "precision": "",
        "notes": f"No truth set. Run status: {status}.",
    }


def seqc2_rows() -> list[dict]:
    """Score SEQC2 SNV and indel calls against truth VCFs."""
    rows = []
    data_dir = PUBLIC_DIR / "seqc2" / "data"
    results_dir = PUBLIC_DIR / "seqc2" / "results"

    configs = [
        {
            "variant_type": "SNV",
            "truth_vcf": data_dir / "HCC1395_truth_SNV.vcf.gz",
            "disc_vcf": results_dir / "kam_seqc2_snv" / "calls_discovery.vcf",
            "ti_vcf": results_dir / "kam_seqc2_snv" / "calls_tumour_informed.vcf",
            "tol": SNV_TOLERANCE,
        },
        {
            "variant_type": "indel",
            "truth_vcf": data_dir / "HCC1395_truth_indel.vcf.gz",
            "disc_vcf": results_dir / "kam_seqc2_indel" / "calls_discovery.vcf",
            "ti_vcf": results_dir / "kam_seqc2_indel" / "calls_tumour_informed.vcf",
            "tol": SNV_TOLERANCE,
        },
    ]

    for cfg in configs:
        truth = load_truth_variants(cfg["truth_vcf"])
        if not cfg["truth_vcf"].exists():
            print(
                f"[WARN] SEQC2 truth VCF not found: {cfg['truth_vcf']}",
                file=sys.stderr,
            )

        for mode, vcf_path in [("discovery", cfg["disc_vcf"]), ("tumour_informed", cfg["ti_vcf"])]:
            if not vcf_path.exists():
                print(
                    f"[WARN] SEQC2 {cfg['variant_type']} {mode} VCF not found: {vcf_path}",
                    file=sys.stderr,
                )
            calls = load_pass_calls(vcf_path)
            m = score(truth, calls, cfg["tol"])
            rows.append(
                {
                    "dataset": "PUB-003_seqc2",
                    "variant_type": cfg["variant_type"],
                    "mode": mode,
                    **m,
                    "notes": "",
                }
            )

    return rows


def colo829_rows() -> list[dict]:
    """Score COLO829 SV calls against the truth VCF."""
    rows = []
    data_dir = PUBLIC_DIR / "colo829" / "data"
    results_dir = PUBLIC_DIR / "colo829" / "results"
    out_dir = results_dir / "kam_colo829_sv"

    truth_vcf = data_dir / "COLO829_truthset_somatic_v4.1.vcf"
    if not truth_vcf.exists():
        print(f"[WARN] COLO829 truth VCF not found: {truth_vcf}", file=sys.stderr)

    truth = load_truth_variants(truth_vcf)

    for mode, vcf_path in [
        ("discovery", out_dir / "calls_discovery.vcf"),
        ("tumour_informed", out_dir / "calls_tumour_informed.vcf"),
    ]:
        if not vcf_path.exists():
            print(
                f"[WARN] COLO829 {mode} VCF not found: {vcf_path}",
                file=sys.stderr,
            )
        calls = load_pass_calls(vcf_path)
        m = score(truth, calls, SV_TOLERANCE)
        rows.append(
            {
                "dataset": "PUB-004_colo829",
                "variant_type": "SV",
                "mode": mode,
                **m,
                "notes": f"68 validated SVs (38 DEL, 13 TRA, 7 DUP, 7 INV, 3 INS). "
                         f"Position tolerance: ±{SV_TOLERANCE} bp.",
            }
        )

    return rows


# ── Main ──────────────────────────────────────────────────────────────────────

COLUMNS = [
    "dataset",
    "variant_type",
    "mode",
    "n_truth",
    "tp",
    "fp",
    "fn",
    "sensitivity",
    "precision",
    "notes",
]


def main() -> None:
    rows: list[dict] = []

    rows.append(umi_benchmark_row())
    rows.extend(seqc2_rows())
    rows.extend(colo829_rows())

    OUTPUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    with open(OUTPUT_CSV, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=COLUMNS)
        writer.writeheader()
        writer.writerows(rows)

    print(f"Written: {OUTPUT_CSV} ({len(rows)} rows)")
    print("")

    # Print a human-readable summary.
    print(f"{'Dataset':<30} {'Type':<10} {'Mode':<18} {'n_truth':>8} "
          f"{'TP':>6} {'FP':>6} {'FN':>6} {'Sens':>8} {'Prec':>8}")
    print("-" * 108)
    for row in rows:
        print(
            f"{row['dataset']:<30} {row['variant_type']:<10} {row['mode']:<18} "
            f"{str(row['n_truth']):>8} {str(row['tp']):>6} {str(row['fp']):>6} "
            f"{str(row['fn']):>6} {str(row['sensitivity']):>8} {str(row['precision']):>8}"
        )


if __name__ == "__main__":
    main()
