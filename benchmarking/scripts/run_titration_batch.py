#!/usr/bin/env python3
"""Run kam on all titration samples and score against truth variants.

Uses up to READS_PER_SAMPLE read pairs per sample. A per-sample RSS cap
(--rss-limit-gb) kills any sample that exceeds the budget mid-run and
marks it as skipped rather than crashing the whole batch.
Anonymises all sample identifiers in output.
"""

import argparse
import os
import re
import subprocess
import threading
import queue
import time
import tempfile
import shutil
import csv
import sys
import psutil
from pathlib import Path

# ── Defaults ──────────────────────────────────────────────────────────────────
REPO = Path(__file__).resolve().parents[2]
_DEFAULT_KAM     = REPO / "target/release/kam"
_DEFAULT_TARGETS = REPO / "benchmarking/scripts/targets_100bp.fa"
_DEFAULT_TRUTH   = REPO / "benchmarking/scripts/truth_variants.vcf"
_DEFAULT_FASTQ   = Path(os.environ.get("KAM_FASTQ_DIR", "/data/titration-nondedup/fastqs"))
_DEFAULT_RESULTS = REPO / "benchmarking/results/tables"

# These are overridden by parse_args() at startup.
KAM          = _DEFAULT_KAM
TARGETS      = _DEFAULT_TARGETS
TRUTH_VCF    = _DEFAULT_TRUTH
FASTQ_DIR    = _DEFAULT_FASTQ
RESULTS_DIR  = _DEFAULT_RESULTS
RESULTS_FILE = RESULTS_DIR / "titration_results.tsv"

READS_PER_SAMPLE = 1_000_000  # 1M read pairs; most samples peak ~1.6 GB
PEAK_RSS_LIMIT_MB = 5 * 1024  # kill and skip any sample that exceeds 5 GB
TRUTH_PANEL_SIZE = 375        # total truth variants in the panel

# Optional caller flags (set to non-None to pass them to kam run).
KAM_MAX_VAF: float | None = None
KAM_MIN_ALT_MOLECULES: int | None = None
KAM_MIN_FAMILY_SIZE: int | None = None
KAM_TARGET_VARIANTS: Path | None = None

# ── Truth variants ─────────────────────────────────────────────────────────────
def load_truth_set(vcf_path):
    """Load truth variants as (chrom, pos, ref, alt) tuples for positional matching.

    Returns:
        truth_all:   full set of all truth variants
        truth_snv:   SNP/MNV variants only
        truth_indel: insertion/deletion variants only
    """
    truth_all = set()
    truth_snv = set()
    truth_indel = set()
    with open(vcf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 8:
                continue
            chrom, pos, _, ref, alt, _, _, info = (
                parts[0], int(parts[1]), parts[2], parts[3], parts[4],
                parts[5], parts[6], parts[7],
            )
            key = (chrom, pos, ref, alt)
            truth_all.add(key)
            vtype = "SNP" if "TYPE=SNP" in info else "INDEL"
            if vtype == "SNP":
                truth_snv.add(key)
            else:
                truth_indel.add(key)
    return truth_all, truth_snv, truth_indel


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
            # Trim common suffix to get the minimal differing region.
            common_suffix = 0
            while (common_suffix < len(ref_seq) - diff_pos - 1
                   and common_suffix < len(alt_seq) - diff_pos - 1
                   and ref_seq[-(common_suffix + 1)] == alt_seq[-(common_suffix + 1)]):
                common_suffix += 1

            if common_suffix > 0:
                ref_trimmed = ref_seq[diff_pos:-common_suffix]
                alt_trimmed = alt_seq[diff_pos:-common_suffix]
            else:
                ref_trimmed = ref_seq[diff_pos:]
                alt_trimmed = alt_seq[diff_pos:]

            if len(ref_seq) == len(alt_seq):
                # SNV: single-base substitution.
                ref_allele = ref_trimmed[0] if len(ref_trimmed) == 1 else ref_trimmed
                alt_allele = alt_trimmed[0] if len(alt_trimmed) == 1 else alt_trimmed
                genomic_pos = target_start + diff_pos
            else:
                # Indel: find the minimal allele pair, then left-normalise to
                # match the VCF convention used in the truth set.

                # Remove any inner common suffix between ref_trimmed and
                # alt_trimmed (the outer common suffix was already stripped).
                # Use len() without -1 so alt_min can reach empty (pure deletion).
                inner_cs = 0
                while (inner_cs < len(ref_trimmed)
                       and inner_cs < len(alt_trimmed)
                       and ref_trimmed[-(inner_cs + 1)] == alt_trimmed[-(inner_cs + 1)]):
                    inner_cs += 1
                ref_min = ref_trimmed[:-inner_cs] if inner_cs > 0 else ref_trimmed
                alt_min = alt_trimmed[:-inner_cs] if inner_cs > 0 else alt_trimmed

                # Remove any inner common prefix as well.
                inner_cp = 0
                for r, a in zip(ref_min, alt_min):
                    if r != a:
                        break
                    inner_cp += 1
                ref_min = ref_min[inner_cp:]
                alt_min = alt_min[inner_cp:]
                indel_start = diff_pos + inner_cp  # index of first ins/del base

                if len(ref_seq) > len(alt_seq):
                    # Deletion: ref_min is the deleted sequence.
                    del_seq = ref_min
                    anchor_pos = indel_start - 1

                    # Left-normalise: shift anchor left while the preceding
                    # base matches the last base of the deleted sequence.
                    # This converts right-aligned deletions in repeat runs to
                    # the left-aligned representation used by VCF.
                    while (anchor_pos > 0
                           and del_seq
                           and ref_seq[anchor_pos] == del_seq[-1]):
                        del_seq = del_seq[-1:] + del_seq[:-1]  # rotate right
                        anchor_pos -= 1

                    if anchor_pos >= 0:
                        anchor = ref_seq[anchor_pos]
                        ref_allele = anchor + del_seq
                        alt_allele = anchor
                        genomic_pos = target_start + anchor_pos
                    else:
                        ref_allele = del_seq
                        alt_allele = ""
                        genomic_pos = target_start + indel_start
                else:
                    # Insertion: alt_min is the inserted sequence.
                    ins_seq = alt_min
                    anchor_pos = indel_start - 1

                    # Left-normalise insertions similarly.
                    while (anchor_pos > 0
                           and ins_seq
                           and ref_seq[anchor_pos] == ins_seq[-1]):
                        ins_seq = ins_seq[-1:] + ins_seq[:-1]  # rotate right
                        anchor_pos -= 1

                    if anchor_pos >= 0:
                        anchor = ref_seq[anchor_pos]
                        ref_allele = anchor
                        alt_allele = anchor + ins_seq
                        genomic_pos = target_start + anchor_pos
                    else:
                        ref_allele = ""
                        alt_allele = ins_seq
                        genomic_pos = target_start + indel_start
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
def _confusion(called, truth_subset, panel_size):
    """Compute confusion-matrix metrics for one (called, truth) pair."""
    tp = len(called & truth_subset)
    fp = len(called - truth_subset)
    fn = len(truth_subset - called)
    tn = max(0, panel_size - tp - fp - fn)
    sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    precision   = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    specificity = tn / (tn + fp) if (tn + fp) > 0 else 0.0
    fpr         = fp / (fp + tn) if (fp + tn) > 0 else 0.0
    fdr         = fp / (tp + fp) if (tp + fp) > 0 else 0.0
    fnr         = fn / (fn + tp) if (fn + tp) > 0 else 0.0
    f1 = (2 * precision * sensitivity / (precision + sensitivity)
          if (precision + sensitivity) > 0 else 0.0)
    return {
        "tp": tp, "fp": fp, "fn": fn, "tn": tn,
        "sensitivity": sensitivity, "precision": precision,
        "specificity": specificity, "fpr": fpr, "fdr": fdr, "fnr": fnr,
        "f1": f1,
    }


def score_tsv(called_tsv, truth_all, truth_snv, truth_indel):
    """Score called variants; return full confusion-matrix metrics overall and by type.

    For a targeted panel with TRUTH_PANEL_SIZE known variant sites:
      TP  = called AND in truth
      FP  = called AND NOT in truth
      FN  = NOT called AND in truth
      TN  = panel sites with no FP call and no missed truth call
    """
    called = extract_called_variants(called_tsv)
    overall = _confusion(called, truth_all, TRUTH_PANEL_SIZE)
    snv     = _confusion(called, truth_snv, len(truth_snv))
    indel   = _confusion(called, truth_indel, len(truth_indel))
    return {"overall": overall, "snv": snv, "indel": indel, "n_called": len(called)}

# Stage tags in the order kam emits them.
_STAGE_TAGS = [
    ("assemble", "[run/assemble]"),
    ("index",    "[run/index]"),
    ("pathfind", "[run/pathfind]"),
    ("call",     "[run/call]"),
    ("output",   "[run] output"),
]

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
    if KAM_MAX_VAF is not None:
        cmd += ["--max-vaf", str(KAM_MAX_VAF)]
    if KAM_MIN_ALT_MOLECULES is not None:
        cmd += ["--min-alt-molecules", str(KAM_MIN_ALT_MOLECULES)]
    if KAM_MIN_FAMILY_SIZE is not None:
        cmd += ["--min-family-size", str(KAM_MIN_FAMILY_SIZE)]
    if KAM_TARGET_VARIANTS is not None:
        cmd += ["--target-variants", str(KAM_TARGET_VARIANTS)]

    print(f"  [{name}] running kam...", flush=True)
    t0 = time.time()
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, text=True)

    # ── Per-stage RSS/CPU tracking ────────────────────────────────────────────
    # Read stdout in a background thread; RSS is sampled in the main loop.
    # When a stage-completion log line arrives, the sampler records peak RSS
    # and CPU for the stage that just finished.
    stdout_q: queue.Queue[str | None] = queue.Queue()

    def _reader():
        for ln in proc.stdout:
            stdout_q.put(ln)
        stdout_q.put(None)

    threading.Thread(target=_reader, daemon=True).start()

    # Per-stage accumulators.
    stage_peak_rss  = {s: 0.0 for s, _ in _STAGE_TAGS}
    stage_peak_cpu  = {s: 0.0 for s, _ in _STAGE_TAGS}
    # Stages complete in order; current_stage is the one we are in now.
    stage_names = [s for s, _ in _STAGE_TAGS]
    stage_idx   = 0  # index into stage_names; current stage = stage_names[stage_idx]

    stdout_lines: list[str] = []
    peak_mb = 0.0
    killed_oom = False
    ps_proc = psutil.Process(proc.pid)

    while True:
        # Drain available stdout lines, advancing stage pointer when log tags appear.
        try:
            while True:
                ln = stdout_q.get_nowait()
                if ln is None:
                    break
                stdout_lines.append(ln)
                _, tag = _STAGE_TAGS[stage_idx] if stage_idx < len(_STAGE_TAGS) else (None, None)
                if tag and tag in ln:
                    # Stage completed; advance to next stage.
                    stage_idx = min(stage_idx + 1, len(stage_names) - 1)
        except queue.Empty:
            pass

        # Sample RSS and CPU.
        try:
            mem  = ps_proc.memory_info().rss / 1024 / 1024
            cpu  = ps_proc.cpu_percent(interval=None)
            peak_mb = max(peak_mb, mem)
            if stage_idx < len(stage_names):
                s = stage_names[stage_idx]
                stage_peak_rss[s] = max(stage_peak_rss[s], mem)
                stage_peak_cpu[s] = max(stage_peak_cpu[s], cpu)
            if mem > PEAK_RSS_LIMIT_MB:
                proc.kill()
                killed_oom = True
                break
        except Exception:
            pass

        if proc.poll() is not None:
            # Process exited; drain remaining lines.
            try:
                while True:
                    ln = stdout_q.get(timeout=0.2)
                    if ln is None:
                        break
                    stdout_lines.append(ln)
            except queue.Empty:
                pass
            break

        time.sleep(0.1)

    proc.wait()
    elapsed = time.time() - t0

    if killed_oom:
        print(f"  [{name}] KILLED: exceeded {PEAK_RSS_LIMIT_MB / 1024:.1f} GB RSS at {peak_mb / 1024:.2f} GB",
              flush=True)
        return None  # caller will write a skipped row

    stdout = "".join(stdout_lines)

    # Parse stage stats and per-stage timing from stdout.
    molecules = duplex = variants_pass = 0
    stage_ms = {"assemble": 0, "index": 0, "pathfind": 0, "call": 0, "output": 0}
    for line in stdout.splitlines():
        if "[run/assemble]" in line:
            for key in ("molecules", "duplex"):
                m = re.search(rf"{key}=(\d+)", line)
                if m:
                    if key == "molecules":
                        molecules = int(m.group(1))
                    else:
                        duplex = int(m.group(1))
            m = re.search(r"time_ms=(\d+)", line)
            if m:
                stage_ms["assemble"] = int(m.group(1))
        elif "[run/index]" in line:
            m = re.search(r"time_ms=(\d+)", line)
            if m:
                stage_ms["index"] = int(m.group(1))
        elif "[run/pathfind]" in line:
            m = re.search(r"time_ms=(\d+)", line)
            if m:
                stage_ms["pathfind"] = int(m.group(1))
        elif "[run/call]" in line:
            m = re.search(r"pass=(\d+)", line)
            if m:
                variants_pass = int(m.group(1))
            m = re.search(r"time_ms=(\d+)", line)
            if m:
                stage_ms["call"] = int(m.group(1))
        elif "[run] output" in line:
            m = re.search(r"time_ms=(\d+)", line)
            if m:
                stage_ms["output"] = int(m.group(1))

    called_tsv = out_dir / "variants.tsv"
    scores = score_tsv(called_tsv, *truth_set)
    ov = scores["overall"]
    sv = scores["snv"]
    ind = scores["indel"]

    def r(v): return round(v, 4)

    row = {
        "sample": name,
        "ng": sample["ng"],
        "vaf": sample["vaf"],
        "molecules": molecules,
        "duplex": duplex,
        "duplex_pct": round(100 * duplex / molecules, 3) if molecules else 0,
        "variants_called": variants_pass,
        # overall
        "tp": ov["tp"], "fp": ov["fp"], "fn": ov["fn"], "tn": ov["tn"],
        "sensitivity": r(ov["sensitivity"]), "precision": r(ov["precision"]),
        "specificity": r(ov["specificity"]), "fpr": r(ov["fpr"]),
        "fdr": r(ov["fdr"]), "fnr": r(ov["fnr"]), "f1": r(ov["f1"]),
        # SNV subset
        "snv_tp": sv["tp"], "snv_fp": sv["fp"], "snv_fn": sv["fn"],
        "snv_sensitivity": r(sv["sensitivity"]), "snv_precision": r(sv["precision"]),
        "snv_specificity": r(sv["specificity"]), "snv_fpr": r(sv["fpr"]),
        # indel subset
        "indel_tp": ind["tp"], "indel_fp": ind["fp"], "indel_fn": ind["fn"],
        "indel_sensitivity": r(ind["sensitivity"]), "indel_precision": r(ind["precision"]),
        "indel_specificity": r(ind["specificity"]), "indel_fpr": r(ind["fpr"]),
        # per-stage timing (ms)
        "t_assemble_ms": stage_ms["assemble"],
        "t_index_ms": stage_ms["index"],
        "t_pathfind_ms": stage_ms["pathfind"],
        "t_call_ms": stage_ms["call"],
        "t_output_ms": stage_ms["output"],
        # per-stage peak RSS (MB) and peak CPU (%)
        "rss_assemble_mb": round(stage_peak_rss["assemble"], 1),
        "rss_index_mb":    round(stage_peak_rss["index"],    1),
        "rss_pathfind_mb": round(stage_peak_rss["pathfind"], 1),
        "rss_call_mb":     round(stage_peak_rss["call"],     1),
        "rss_output_mb":   round(stage_peak_rss["output"],   1),
        "cpu_assemble_pct": round(stage_peak_cpu["assemble"], 1),
        "cpu_index_pct":    round(stage_peak_cpu["index"],    1),
        "cpu_pathfind_pct": round(stage_peak_cpu["pathfind"], 1),
        "cpu_call_pct":     round(stage_peak_cpu["call"],     1),
        "cpu_output_pct":   round(stage_peak_cpu["output"],   1),
        # runtime
        "peak_rss_mb": round(peak_mb, 1),
        "wall_time_s": round(elapsed, 1),
        "exit_code": proc.returncode,
    }

    status = "OK" if proc.returncode == 0 else "FAIL"
    print(f"  [{name}] {status} | "
          f"mols={molecules:,} | "
          f"sens={ov['sensitivity']:.3f} spec={ov['specificity']:.3f} "
          f"prec={ov['precision']:.3f} | "
          f"snv={sv['sensitivity']:.3f} indel={ind['sensitivity']:.3f} | "
          f"assemble={stage_ms['assemble']}ms({stage_peak_rss['assemble']:.0f}MB) "
          f"index={stage_ms['index']}ms({stage_peak_rss['index']:.0f}MB) "
          f"pathfind={stage_ms['pathfind']}ms({stage_peak_rss['pathfind']:.0f}MB) "
          f"call={stage_ms['call']}ms({stage_peak_rss['call']:.0f}MB) | "
          f"{elapsed:.1f}s peak={peak_mb:.0f}MB", flush=True)
    return row

# ── Main ──────────────────────────────────────────────────────────────────────
def main():
    global KAM, TARGETS, TRUTH_VCF, FASTQ_DIR, RESULTS_DIR, RESULTS_FILE
    global READS_PER_SAMPLE, PEAK_RSS_LIMIT_MB
    global KAM_MAX_VAF, KAM_MIN_ALT_MOLECULES, KAM_MIN_FAMILY_SIZE, KAM_TARGET_VARIANTS

    parser = argparse.ArgumentParser(
        description="Run kam on all titration samples and score against truth variants."
    )
    parser.add_argument("--fastq-dir", type=Path, default=_DEFAULT_FASTQ,
                        help="Directory containing TWIST_STDV2_* FASTQ pairs "
                             "(default: $KAM_FASTQ_DIR or /data/titration-nondedup/fastqs)")
    parser.add_argument("--truth-vcf", type=Path, default=_DEFAULT_TRUTH,
                        help="Truth VCF file (default: benchmarking/scripts/truth_variants.vcf)")
    parser.add_argument("--targets", type=Path, default=_DEFAULT_TARGETS,
                        help="Target FASTA file (default: benchmarking/scripts/targets_100bp.fa)")
    parser.add_argument("--kam-binary", type=Path, default=_DEFAULT_KAM,
                        help="Path to the kam binary (default: target/release/kam)")
    parser.add_argument("--results-dir", type=Path, default=_DEFAULT_RESULTS,
                        help="Directory for output TSV (default: benchmarking/results/tables)")
    parser.add_argument("--output", type=str, default=None,
                        help="Output TSV filename within --results-dir "
                             "(default: titration_results_{N}reads.tsv based on --reads)")
    parser.add_argument("--reads", type=int, default=READS_PER_SAMPLE,
                        help=f"Read pairs per sample (default: {READS_PER_SAMPLE:,})")
    parser.add_argument("--rss-limit-gb", type=float, default=PEAK_RSS_LIMIT_MB / 1024,
                        help=f"Peak RSS kill threshold in GB (default: {PEAK_RSS_LIMIT_MB / 1024:.0f})")
    parser.add_argument("--max-vaf", type=float, default=None,
                        help="Maximum VAF for a PASS call; calls above this are labelled HighVaf. "
                             "Use 0.35 to exclude germline heterozygous variants.")
    parser.add_argument("--min-alt-molecules", type=int, default=None,
                        help="Minimum alt-supporting molecules to emit a call (default: 2).")
    parser.add_argument("--min-family-size", type=int, default=None,
                        help="Minimum reads per UMI family to keep a molecule (default: 1). "
                             "Use 2 to replicate HUMID-style singleton filtering.")
    parser.add_argument("--target-variants", type=Path, default=None,
                        help="VCF of expected somatic variants for tumour-informed monitoring "
                             "mode. Only calls matching (CHROM, POS, REF, ALT) in this VCF "
                             "are marked PASS. Suppresses background biological FPs.")
    args = parser.parse_args()

    KAM             = args.kam_binary
    TARGETS         = args.targets
    TRUTH_VCF       = args.truth_vcf
    FASTQ_DIR       = args.fastq_dir
    RESULTS_DIR     = args.results_dir
    READS_PER_SAMPLE  = args.reads
    PEAK_RSS_LIMIT_MB = int(args.rss_limit_gb * 1024)
    KAM_MAX_VAF           = args.max_vaf
    KAM_MIN_ALT_MOLECULES = args.min_alt_molecules
    KAM_MIN_FAMILY_SIZE   = args.min_family_size
    KAM_TARGET_VARIANTS   = args.target_variants

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    if args.output:
        output_name = args.output
    else:
        n = READS_PER_SAMPLE
        if n % 1_000_000 == 0:
            label = f"{n // 1_000_000}m"
        else:
            label = f"{n // 1_000}k"
        output_name = f"titration_results_{label}reads.tsv"
    RESULTS_FILE = RESULTS_DIR / output_name

    truth_set = load_truth_set(TRUTH_VCF)
    truth_all, truth_snv, truth_indel = truth_set
    print(f"Truth variants loaded: {len(truth_all)} total "
          f"({len(truth_snv)} SNV, {len(truth_indel)} indel)", flush=True)

    samples = find_samples()
    print(f"Found {len(samples)} samples", flush=True)

    fieldnames = [
        "sample", "ng", "vaf", "molecules", "duplex", "duplex_pct",
        "variants_called",
        # overall confusion matrix
        "tp", "fp", "fn", "tn",
        "sensitivity", "precision", "specificity", "fpr", "fdr", "fnr", "f1",
        # SNV subset
        "snv_tp", "snv_fp", "snv_fn",
        "snv_sensitivity", "snv_precision", "snv_specificity", "snv_fpr",
        # indel subset
        "indel_tp", "indel_fp", "indel_fn",
        "indel_sensitivity", "indel_precision", "indel_specificity", "indel_fpr",
        # per-stage timing (ms)
        "t_assemble_ms", "t_index_ms", "t_pathfind_ms", "t_call_ms", "t_output_ms",
        # per-stage peak RSS (MB) and peak CPU (%)
        "rss_assemble_mb", "rss_index_mb", "rss_pathfind_mb", "rss_call_mb", "rss_output_mb",
        "cpu_assemble_pct", "cpu_index_pct", "cpu_pathfind_pct", "cpu_call_pct", "cpu_output_pct",
        # overall runtime
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
                    row = None
                if row is None:
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
