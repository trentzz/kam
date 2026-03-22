#!/usr/bin/env python3
"""Generate per-sample discovery and monitoring TSV files for all titration samples.

For each sample, runs kam twice:
  1. Discovery mode: no --target-variants. All quality-passing calls included.
  2. Monitoring mode: --target-variants truth_variants.vcf. Only calls matching
     the truth somatic panel are labelled PASS; all others are NotTargeted.

Output: two TSV files per sample in docs/benchmarking/snvindel/per_sample/:
  <sample>_discovery.tsv
  <sample>_monitoring.tsv

Timing for every run is appended to timing.tsv in the same directory:
  sample  mode  n_reads  total_ms  assemble_ms  index_ms  pathfind_ms  call_ms

The script skips any sample/mode whose TSV AND timing row both already exist,
so it is safe to rerun after partial completion.

Requires KAM_FASTQ_DIR to point to the titration FASTQ directory.

Usage:
    python generate_per_sample_tsvs.py [--reads N] [--output-dir DIR]
"""

import argparse
import csv
import os
import re
import subprocess
import sys
import tempfile
import shutil
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
KAM          = REPO / "target/release/kam"
TARGETS      = REPO / "docs/benchmarking/scripts/targets_100bp.fa"
TRUTH_VCF    = REPO / "docs/benchmarking/scripts/truth_variants.vcf"
_DEFAULT_OUT = REPO / "docs/benchmarking/snvindel/per_sample"
_DEFAULT_FASTQ = Path(os.environ.get("KAM_FASTQ_DIR", "/data/titration-nondedup/fastqs"))

PATTERN = re.compile(r"TWIST_STDV2_(\d+ng)_VAF_(\w+)pc_.*_R1\.fastq\.gz$")
READS_PER_SAMPLE = 2_000_000

TIMING_COLS = ["sample", "mode", "n_reads", "total_ms",
               "assemble_ms", "index_ms", "pathfind_ms", "call_ms"]

# Regexes for stage timings in kam stderr.
_RE = {
    "assemble_ms":  re.compile(r"\[run/assemble\].*time_ms=(\d+)"),
    "index_ms":     re.compile(r"\[run/index\].*time_ms=(\d+)"),
    "pathfind_ms":  re.compile(r"\[run/pathfind\].*time_ms=(\d+)"),
    "call_ms":      re.compile(r"\[run/call\].*time_ms=(\d+)"),
    "total_ms":     re.compile(r"\[run\] output.*total_ms=(\d+)"),
}


def vaf_label_to_float(label: str) -> float:
    return float(label.replace("p", "."))


def find_samples(fastq_dir: Path) -> list[dict]:
    samples = []
    for f in sorted(fastq_dir.glob("*_R1.fastq.gz")):
        m = PATTERN.match(f.name)
        if not m:
            continue
        ng = m.group(1)
        vaf_label = m.group(2)
        r2 = Path(str(f).replace("_R1.fastq.gz", "_R2.fastq.gz"))
        if not r2.exists():
            print(f"WARNING: no R2 for {f.name}", file=sys.stderr)
            continue
        samples.append({
            "r1": f, "r2": r2,
            "ng": ng, "vaf_label": vaf_label,
            "vaf": vaf_label_to_float(vaf_label),
            "name": f"Sample_{ng}_VAF_{vaf_label}pc",
        })
    return samples


def load_timing_index(timing_path: Path) -> set[tuple[str, str]]:
    """Return set of (sample, mode) pairs already recorded in timing.tsv."""
    if not timing_path.exists():
        return set()
    with open(timing_path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        return {(row["sample"], row["mode"]) for row in reader}


def append_timing(timing_path: Path, row: dict) -> None:
    write_header = not timing_path.exists()
    with open(timing_path, "a") as f:
        writer = csv.DictWriter(f, fieldnames=TIMING_COLS, delimiter="\t")
        if write_header:
            writer.writeheader()
        writer.writerow(row)


def parse_timing(stderr: str) -> dict:
    result = {}
    for key, pat in _RE.items():
        m = pat.search(stderr)
        result[key] = int(m.group(1)) if m else ""
    return result


def subset_fastq(gz_path: Path, n_reads: int, out_path: Path) -> None:
    n_lines = n_reads * 4
    with open(out_path, "wb") as f_out:
        zcat = subprocess.Popen(["zcat", str(gz_path)], stdout=subprocess.PIPE)
        count = 0
        for line in zcat.stdout:
            if count >= n_lines:
                break
            f_out.write(line)
            count += 1
        zcat.stdout.close()
        zcat.wait()


def run_kam(r1: Path, r2: Path, out_dir: Path,
            target_variants: Path | None) -> tuple[Path | None, str]:
    """Run kam. Returns (path to variants.tsv or None, stderr text)."""
    out_dir.mkdir(parents=True, exist_ok=True)
    cmd = [
        str(KAM), "run",
        "--r1", str(r1),
        "--r2", str(r2),
        "--targets", str(TARGETS),
        "--output-dir", str(out_dir),
        "--output-format", "tsv",
    ]
    if target_variants is not None:
        cmd += ["--target-variants", str(target_variants)]

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"  ERROR: kam exited {result.returncode}", file=sys.stderr)
        print(result.stderr[-1000:], file=sys.stderr)
        return None, result.stderr

    tsv = out_dir / "variants.tsv"
    return (tsv if tsv.exists() else None), result.stderr


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate per-sample discovery and monitoring TSVs")
    parser.add_argument("--reads", type=int, default=READS_PER_SAMPLE,
                        help=f"Read pairs per sample (default: {READS_PER_SAMPLE:,})")
    parser.add_argument("--output-dir", type=Path, default=_DEFAULT_OUT,
                        help="Output directory for TSV files")
    parser.add_argument("--fastq-dir", type=Path, default=_DEFAULT_FASTQ,
                        help="Directory containing titration FASTQ files")
    args = parser.parse_args()

    if not KAM.exists():
        print(f"ERROR: kam binary not found at {KAM}. Run: cargo build --release", file=sys.stderr)
        sys.exit(1)
    if not args.fastq_dir.exists():
        print(f"ERROR: FASTQ dir not found: {args.fastq_dir}", file=sys.stderr)
        print("Set KAM_FASTQ_DIR or pass --fastq-dir", file=sys.stderr)
        sys.exit(1)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    timing_path = args.output_dir / "timing.tsv"
    timing_done = load_timing_index(timing_path)

    samples = find_samples(args.fastq_dir)
    if not samples:
        print(f"ERROR: no samples found in {args.fastq_dir}", file=sys.stderr)
        sys.exit(1)

    print(f"Found {len(samples)} samples. Generating discovery + monitoring TSVs...\n")

    for i, sample in enumerate(samples, 1):
        name = sample["name"]
        print(f"[{i}/{len(samples)}] {name}")

        discovery_out = args.output_dir / f"{name}_discovery.tsv"
        monitoring_out = args.output_dir / f"{name}_monitoring.tsv"

        need_disc = not discovery_out.exists() or (name, "discovery") not in timing_done
        need_mon  = not monitoring_out.exists() or (name, "monitoring") not in timing_done

        if not need_disc and not need_mon:
            print(f"  already complete, skipping")
            continue

        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            r1_sub = tmp / "r1.fastq"
            r2_sub = tmp / "r2.fastq"

            print(f"  subsetting {args.reads:,} reads...")
            subset_fastq(sample["r1"], args.reads, r1_sub)
            subset_fastq(sample["r2"], args.reads, r2_sub)

            # Discovery run
            if need_disc:
                print(f"  running discovery mode...")
                disc_dir = tmp / "discovery"
                tsv, stderr = run_kam(r1_sub, r2_sub, disc_dir, target_variants=None)
                if tsv:
                    shutil.copy(tsv, discovery_out)
                    timing = parse_timing(stderr)
                    append_timing(timing_path, {
                        "sample": name, "mode": "discovery",
                        "n_reads": args.reads, **timing,
                    })
                    print(f"  -> {discovery_out.name}  ({timing.get('total_ms', '?')} ms)")
                else:
                    print(f"  FAILED: discovery mode")

            # Monitoring run
            if need_mon:
                print(f"  running monitoring mode...")
                mon_dir = tmp / "monitoring"
                tsv, stderr = run_kam(r1_sub, r2_sub, mon_dir, target_variants=TRUTH_VCF)
                if tsv:
                    shutil.copy(tsv, monitoring_out)
                    timing = parse_timing(stderr)
                    append_timing(timing_path, {
                        "sample": name, "mode": "monitoring",
                        "n_reads": args.reads, **timing,
                    })
                    print(f"  -> {monitoring_out.name}  ({timing.get('total_ms', '?')} ms)")
                else:
                    print(f"  FAILED: monitoring mode")

        print()

    print(f"Done. TSVs written to: {args.output_dir}")
    disc_count = len(list(args.output_dir.glob("*_discovery.tsv")))
    mon_count  = len(list(args.output_dir.glob("*_monitoring.tsv")))
    print(f"  {disc_count} discovery TSVs, {mon_count} monitoring TSVs")
    if timing_path.exists():
        timing_rows = load_timing_index(timing_path)
        print(f"  {len(timing_rows)} timing rows in {timing_path.name}")


if __name__ == "__main__":
    main()
