#!/usr/bin/env python3
"""Run kam on all titration samples and score against truth variants.

Uses 1M read subsets for memory safety (peak ~1.6 GB per sample).
Anonymises all sample identifiers in output.
"""

import os
import re
import subprocess
import time
import tempfile
import shutil
import csv
import sys
import psutil
from pathlib import Path

# ── Paths ─────────────────────────────────────────────────────────────────────
REPO = Path(__file__).resolve().parents[2]
KAM = REPO / "target/release/kam"
TARGETS = REPO / "benchmarking/scripts/targets_100bp.fa"
TRUTH_VCF = REPO / "benchmarking/scripts/truth_variants.vcf"
FASTQ_DIR = Path("/mnt/tzeng-local/tzeng-thesis/titration-nondedup/fastqs")
RESULTS_DIR = REPO / "benchmarking/results/tables"
SCORE_SCRIPT = REPO / "benchmarking/scripts/score_variants.py"

READS_PER_SAMPLE = 200_000  # 200K read pairs → ~350–500 MB peak RSS

RESULTS_DIR.mkdir(parents=True, exist_ok=True)
RESULTS_FILE = RESULTS_DIR / "titration_results.tsv"

# ── Truth variants ─────────────────────────────────────────────────────────────
def load_truth_set(vcf_path):
    """Load truth variants as (chrom, pos, ref, alt) tuples for positional matching."""
    truth = set()
    with open(vcf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) >= 5:
                chrom, pos, _, ref, alt = parts[0], int(parts[1]), parts[2], parts[3], parts[4]
                truth.add((chrom, pos, ref, alt))
    return truth


def extract_called_variants(tsv_path):
    """Parse kam TSV output; derive (chrom, pos, ref, alt) from target_id and sequences.

    target_id is 'chrN:start-end' (0-based start, end exclusive, matching BED).
    ref_seq and alt_seq are the full 100bp sequences; the variant is where they
    first differ (for SNVs) or where a gap/insertion is introduced (for indels).
    We return only PASS variants.
    """
    called = set()
    if not tsv_path.exists():
        return called
    with open(tsv_path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            if row.get("filter") != "PASS":
                continue
            tid = row["target_id"]
            # Parse target_id: chrom:start-end (start is 0-based)
            m = re.match(r"^(chr\w+):(\d+)-(\d+)$", tid)
            if not m:
                continue
            chrom = m.group(1)
            target_start = int(m.group(2))  # 0-based
            ref_seq = row["ref_seq"]
            alt_seq = row["alt_seq"]
            # Find the leftmost differing position
            diff_pos = next(
                (i for i, (r, a) in enumerate(zip(ref_seq, alt_seq)) if r != a),
                None,
            )
            if diff_pos is None:
                continue
            # Determine REF and ALT alleles from the diff
            if len(ref_seq) == len(alt_seq):
                # SNV or MNV: same-length substitution
                ref_allele = ref_seq[diff_pos]
                alt_allele = alt_seq[diff_pos]
                # FASTA headers use 1-based start; diff_pos is 0-based offset → add directly
                genomic_pos = target_start + diff_pos
            else:
                # Indel: extract alleles from first diff position onward
                ref_allele = ref_seq[diff_pos:]
                alt_allele = alt_seq[diff_pos:]
                genomic_pos = target_start + diff_pos
            called.add((chrom, genomic_pos, ref_allele, alt_allele))
    return called

# ── Sample discovery ──────────────────────────────────────────────────────────
PATTERN = re.compile(
    r"TWIST_STDV2_(\d+ng)_VAF_(\w+)pc_.*_R1\.fastq\.gz$"
)

def vaf_label_to_float(label):
    # "0p001" -> 0.001, "0" -> 0.0, "1" -> 1.0
    return float(label.replace("p", "."))

def find_samples():
    samples = []
    for f in sorted(FASTQ_DIR.glob("*_R1.fastq.gz")):
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
            "ng": ng,
            "vaf_label": vaf_label,
            "vaf": vaf_label_to_float(vaf_label),
            "name": f"Sample_{ng}_VAF_{vaf_label}pc",
        })
    return samples

# ── Subset FASTQ ─────────────────────────────────────────────────────────────
def subset_fastq(gz_path, n_reads, out_path):
    """Stream-decompress and write first n_reads read pairs without loading into memory."""
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

# ── Score variants ────────────────────────────────────────────────────────────
def score_tsv(called_tsv, truth_set):
    """Score called variants against truth using positional (chrom, pos, ref, alt) matching."""
    called = extract_called_variants(called_tsv)
    tp = len(called & truth_set)
    fp = len(called - truth_set)
    fn = len(truth_set - called)
    sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    f1 = (2 * precision * sensitivity / (precision + sensitivity)
          if (precision + sensitivity) > 0 else 0.0)
    return {"tp": tp, "fp": fp, "fn": fn,
            "sensitivity": sensitivity, "precision": precision, "f1": f1,
            "n_called": len(called)}

# ── Run one sample ─────────────────────────────────────────────────────────────
def run_sample(sample, truth_set, tmp_dir):
    name = sample["name"]
    print(f"  [{name}] subsetting {READS_PER_SAMPLE:,} reads...", flush=True)

    r1_sub = tmp_dir / "r1.fastq"
    r2_sub = tmp_dir / "r2.fastq"
    subset_fastq(sample["r1"], READS_PER_SAMPLE, r1_sub)
    subset_fastq(sample["r2"], READS_PER_SAMPLE, r2_sub)

    out_dir = tmp_dir / "kam_out"
    out_dir.mkdir()

    cmd = [
        str(KAM), "run",
        "--r1", str(r1_sub), "--r2", str(r2_sub),
        "--targets", str(TARGETS),
        "--output-dir", str(out_dir),
        "--output-format", "vcf,tsv",
    ]

    print(f"  [{name}] running kam...", flush=True)
    t0 = time.time()
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, text=True)
    peak_mb = 0.0
    while proc.poll() is None:
        try:
            rss = psutil.Process(proc.pid).memory_info().rss / 1024 / 1024
            peak_mb = max(peak_mb, rss)
        except Exception:
            pass
        time.sleep(0.1)
    stdout, _ = proc.communicate()
    elapsed = time.time() - t0

    # Parse stage stats from stdout
    molecules = duplex = variants_pass = 0
    for line in stdout.splitlines():
        if "molecules=" in line:
            m = re.search(r"molecules=(\d+)", line)
            if m:
                molecules = int(m.group(1))
            m = re.search(r"duplex=(\d+)", line)
            if m:
                duplex = int(m.group(1))
        if "pass=" in line:
            m = re.search(r"pass=(\d+)", line)
            if m:
                variants_pass = int(m.group(1))

    called_tsv = out_dir / "variants.tsv"
    scores = score_tsv(called_tsv, truth_set)

    row = {
        "sample": name,
        "ng": sample["ng"],
        "vaf": sample["vaf"],
        "molecules": molecules,
        "duplex": duplex,
        "duplex_pct": round(100 * duplex / molecules, 3) if molecules else 0,
        "variants_called": variants_pass,
        "tp": scores["tp"],
        "fp": scores["fp"],
        "fn": scores["fn"],
        "sensitivity": round(scores["sensitivity"], 4),
        "precision": round(scores["precision"], 4),
        "f1": round(scores["f1"], 4),
        "peak_rss_mb": round(peak_mb, 1),
        "wall_time_s": round(elapsed, 1),
        "exit_code": proc.returncode,
    }

    status = "OK" if proc.returncode == 0 else "FAIL"
    print(f"  [{name}] {status} | "
          f"mols={molecules:,} duplex={duplex} | "
          f"sens={scores['sensitivity']:.3f} prec={scores['precision']:.3f} | "
          f"{elapsed:.1f}s {peak_mb:.0f}MB", flush=True)
    return row

# ── Main ──────────────────────────────────────────────────────────────────────
def main():
    truth_set = load_truth_set(TRUTH_VCF)
    print(f"Truth variants loaded: {len(truth_set)}", flush=True)

    samples = find_samples()
    print(f"Found {len(samples)} samples", flush=True)

    fieldnames = [
        "sample", "ng", "vaf", "molecules", "duplex", "duplex_pct",
        "variants_called", "tp", "fp", "fn",
        "sensitivity", "precision", "f1",
        "peak_rss_mb", "wall_time_s", "exit_code",
    ]

    with open(RESULTS_FILE, "w", newline="") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()

        for i, sample in enumerate(samples, 1):
            print(f"\n[{i}/{len(samples)}] {sample['name']}", flush=True)
            with tempfile.TemporaryDirectory(prefix="kam_titration_") as tmp:
                tmp_dir = Path(tmp)
                try:
                    row = run_sample(sample, truth_set, tmp_dir)
                except Exception as e:
                    print(f"  ERROR: {e}", file=sys.stderr)
                    row = {k: "" for k in fieldnames}
                    row["sample"] = sample["name"]
                    row["ng"] = sample["ng"]
                    row["vaf"] = sample["vaf"]
                    row["exit_code"] = -1
                writer.writerow(row)
                csvfile.flush()

    print(f"\nResults written to {RESULTS_FILE}", flush=True)

if __name__ == "__main__":
    main()
