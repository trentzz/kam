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
REPO = Path(__file__).resolve().parents[4]
_DEFAULT_KAM     = REPO / "target/release/kam"
_DEFAULT_TARGETS = REPO / "docs/benchmarking/01-snvindel/scripts/targets_100bp.fa"
_DEFAULT_TRUTH   = REPO / "docs/benchmarking/01-snvindel/scripts/truth_variants.vcf"
_DEFAULT_FASTQ   = Path(os.environ.get("KAM_FASTQ_DIR", "/data/titration-nondedup/fastqs"))
_DEFAULT_RESULTS = REPO / "docs/benchmarking/07-snvindel-ml-boost-v1/results"

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
KAM_MIN_CONFIDENCE: float | None = None
KAM_MIN_FAMILY_SIZE: int | None = None
KAM_TARGET_VARIANTS: Path | None = None
KAM_KMER_SIZE: int | None = None
KAM_MIN_ALT_DUPLEX: int | None = None
KAM_ML_MODEL: str | None = None
KAM_TI_RESCUE: bool = False

# When set, save per-sample VCFs to this directory.
VCF_SAVE_DIR: Path | None = None

# When set, write per-sample per-target TSV to this directory.
PER_SAMPLE_DIR: Path | None = None

# ── Truth variants ─────────────────────────────────────────────────────────────
def load_truth_set(vcf_path):
    """Load truth variants as (chrom, pos, ref, alt) tuples for positional matching."""
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


def extract_called_variants_with_ml(tsv_path):
    """Parse kam TSV; return two sets: standard PASS calls and ML-filtered calls.

    Returns:
        called_standard: set of (chrom, pos, ref, alt) tuples for PASS calls
        called_ml: set of (chrom, pos, ref, alt) tuples for PASS + ml_filter=ML_PASS calls
        call_details: list of dicts with per-call detail for per-target TSV
    """
    called_standard = set()
    called_ml = set()
    call_details = []

    if not tsv_path.exists():
        return called_standard, called_ml, call_details

    with open(tsv_path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            if row.get("filter") != "PASS":
                continue
            tid = row["target_id"]
            m = re.match(r"^(chr\w+):(\d+)-(\d+)$", tid)
            if not m:
                continue
            chrom = m.group(1)
            target_start = int(m.group(2))
            ref_seq = row["ref_seq"]
            alt_seq = row["alt_seq"]
            diff_pos = next(
                (i for i, (r, a) in enumerate(zip(ref_seq, alt_seq)) if r != a),
                None,
            )
            if diff_pos is None:
                continue
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
                ref_allele = ref_trimmed[0] if len(ref_trimmed) == 1 else ref_trimmed
                alt_allele = alt_trimmed[0] if len(alt_trimmed) == 1 else alt_trimmed
                genomic_pos = target_start + diff_pos
            else:
                inner_cs = 0
                while (inner_cs < len(ref_trimmed) and inner_cs < len(alt_trimmed)
                       and ref_trimmed[-(inner_cs + 1)] == alt_trimmed[-(inner_cs + 1)]):
                    inner_cs += 1
                ref_min = ref_trimmed[:-inner_cs] if inner_cs > 0 else ref_trimmed
                alt_min = alt_trimmed[:-inner_cs] if inner_cs > 0 else alt_trimmed
                inner_cp = 0
                for rv, av in zip(ref_min, alt_min):
                    if rv != av:
                        break
                    inner_cp += 1
                ref_min = ref_min[inner_cp:]
                alt_min = alt_min[inner_cp:]
                indel_start = diff_pos + inner_cp
                if len(ref_seq) > len(alt_seq):
                    del_seq = ref_min
                    anchor_pos = indel_start - 1
                    while anchor_pos > 0 and del_seq and ref_seq[anchor_pos] == del_seq[-1]:
                        del_seq = del_seq[-1:] + del_seq[:-1]
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
                    ins_seq = alt_min
                    anchor_pos = indel_start - 1
                    while anchor_pos > 0 and ins_seq and ref_seq[anchor_pos] == ins_seq[-1]:
                        ins_seq = ins_seq[-1:] + ins_seq[:-1]
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

            key = (chrom, genomic_pos, ref_allele, alt_allele)
            called_standard.add(key)

            ml_filter_val = row.get("ml_filter", ".")
            ml_prob_val = row.get("ml_prob", ".")
            confidence_val = row.get("confidence", ".")
            vaf_val = row.get("vaf", ".")
            if ml_filter_val == "ML_PASS":
                called_ml.add(key)

            call_details.append({
                "key": key,
                "ml_filter": ml_filter_val,
                "ml_prob": ml_prob_val,
                "confidence": confidence_val,
                "called_vaf": vaf_val,
            })

    return called_standard, called_ml, call_details


# ── Sample discovery ──────────────────────────────────────────────────────────
PATTERN = re.compile(
    r"TWIST_STDV2_(\d+ng)_VAF_(\w+)pc_.*_R1\.fastq\.gz$"
)

def vaf_label_to_float(label):
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

    When the TSV contains ml_filter/ml_prob columns (from --ml-model), also computes
    ML-filtered metrics using only calls where ml_filter=ML_PASS.
    """
    called_standard, called_ml, call_details = extract_called_variants_with_ml(called_tsv)
    has_ml = any(d["ml_filter"] != "." for d in call_details)

    overall = _confusion(called_standard, truth_all, TRUTH_PANEL_SIZE)
    snv     = _confusion(called_standard, truth_snv, len(truth_snv))
    indel   = _confusion(called_standard, truth_indel, len(truth_indel))

    if has_ml:
        ml_overall = _confusion(called_ml, truth_all, TRUTH_PANEL_SIZE)
        ml_snv     = _confusion(called_ml, truth_snv, len(truth_snv))
        ml_indel   = _confusion(called_ml, truth_indel, len(truth_indel))
    else:
        ml_overall = ml_snv = ml_indel = None

    return {
        "overall": overall, "snv": snv, "indel": indel,
        "n_called": len(called_standard),
        "ml_overall": ml_overall, "ml_snv": ml_snv, "ml_indel": ml_indel,
        "call_details": call_details,
        "called_standard": called_standard,
        "called_ml": called_ml,
    }


def write_per_target_tsv(sample, scores, truth_all, truth_snv, truth_indel, out_path):
    """Write a per-target TSV with one row per truth variant for this sample.

    Columns: variant_id, chrom, pos, ref, alt, variant_type,
             sample, ng, vaf,
             detected, ml_detected, called_vaf, confidence, ml_prob, ml_filter
    """
    # Build a lookup from variant key to call detail.
    detail_by_key = {d["key"]: d for d in scores["call_details"]}
    called_standard = scores["called_standard"]
    called_ml = scores["called_ml"]
    has_ml = scores["ml_overall"] is not None

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "variant_id", "chrom", "pos", "ref", "alt", "variant_type",
        "sample", "ng", "vaf",
        "detected", "ml_detected", "called_vaf", "confidence", "ml_prob", "ml_filter",
    ]
    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        # Collect and sort truth variants for deterministic output.
        all_truth = sorted(
            [(k, "SNP" if k in truth_snv else "INDEL") for k in truth_all],
            key=lambda x: (x[0][0], x[0][1]),
        )
        for (chrom, pos, ref, alt), vtype in all_truth:
            key = (chrom, pos, ref, alt)
            detected = key in called_standard
            detail = detail_by_key.get(key, {})
            row = {
                "variant_id": f"{chrom}:{pos}:{ref}:{alt}",
                "chrom": chrom,
                "pos": pos,
                "ref": ref,
                "alt": alt,
                "variant_type": vtype,
                "sample": sample["name"],
                "ng": sample["ng"],
                "vaf": sample["vaf"],
                "detected": detected,
                "ml_detected": (key in called_ml) if has_ml else "NA",
                "called_vaf": detail.get("called_vaf", "NA") if detected else "NA",
                "confidence": detail.get("confidence", "NA") if detected else "NA",
                "ml_prob": detail.get("ml_prob", "NA") if detected else "NA",
                "ml_filter": detail.get("ml_filter", "NA") if detected else "NA",
            }
            writer.writerow(row)


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
        "--output-format-override", "vcf,tsv",
    ]
    if KAM_MAX_VAF is not None:
        cmd += ["--max-vaf", str(KAM_MAX_VAF)]
    if KAM_MIN_ALT_MOLECULES is not None:
        cmd += ["--min-alt-molecules", str(KAM_MIN_ALT_MOLECULES)]
    if KAM_MIN_CONFIDENCE is not None:
        cmd += ["--min-confidence", str(KAM_MIN_CONFIDENCE)]
    if KAM_MIN_FAMILY_SIZE is not None:
        cmd += ["--min-family-size", str(KAM_MIN_FAMILY_SIZE)]
    if KAM_TARGET_VARIANTS is not None:
        cmd += ["--target-variants", str(KAM_TARGET_VARIANTS)]
    if KAM_KMER_SIZE is not None:
        cmd += ["-k", str(KAM_KMER_SIZE)]
    if KAM_MIN_ALT_DUPLEX is not None:
        cmd += ["--min-alt-duplex", str(KAM_MIN_ALT_DUPLEX)]
    if KAM_ML_MODEL is not None:
        cmd += ["--ml-model", KAM_ML_MODEL]
    if KAM_TI_RESCUE:
        cmd += ["--ti-rescue"]

    print(f"  [{name}] running kam...", flush=True)
    t0 = time.time()
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, text=True)

    stdout_q: queue.Queue[str | None] = queue.Queue()

    def _reader():
        for ln in proc.stdout:
            stdout_q.put(ln)
        stdout_q.put(None)

    threading.Thread(target=_reader, daemon=True).start()

    stage_peak_rss  = {s: 0.0 for s, _ in _STAGE_TAGS}
    stage_peak_cpu  = {s: 0.0 for s, _ in _STAGE_TAGS}
    stage_names = [s for s, _ in _STAGE_TAGS]
    stage_idx   = 0

    stdout_lines: list[str] = []
    peak_mb = 0.0
    killed_oom = False
    ps_proc = psutil.Process(proc.pid)

    while True:
        try:
            while True:
                ln = stdout_q.get_nowait()
                if ln is None:
                    break
                stdout_lines.append(ln)
                _, tag = _STAGE_TAGS[stage_idx] if stage_idx < len(_STAGE_TAGS) else (None, None)
                if tag and tag in ln:
                    stage_idx = min(stage_idx + 1, len(stage_names) - 1)
        except queue.Empty:
            pass

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
        return None

    stdout = "".join(stdout_lines)

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
    ml_ov  = scores["ml_overall"]
    ml_sv  = scores["ml_snv"]
    ml_ind = scores["ml_indel"]
    has_ml = ml_ov is not None

    def r(v): return round(v, 4)

    row = {
        "sample": name,
        "ng": sample["ng"],
        "vaf": sample["vaf"],
        "molecules": molecules,
        "duplex": duplex,
        "duplex_pct": round(100 * duplex / molecules, 3) if molecules else 0,
        "variants_called": variants_pass,
        "tp": ov["tp"], "fp": ov["fp"], "fn": ov["fn"], "tn": ov["tn"],
        "sensitivity": r(ov["sensitivity"]), "precision": r(ov["precision"]),
        "specificity": r(ov["specificity"]), "fpr": r(ov["fpr"]),
        "fdr": r(ov["fdr"]), "fnr": r(ov["fnr"]), "f1": r(ov["f1"]),
        "snv_tp": sv["tp"], "snv_fp": sv["fp"], "snv_fn": sv["fn"],
        "snv_sensitivity": r(sv["sensitivity"]), "snv_precision": r(sv["precision"]),
        "snv_specificity": r(sv["specificity"]), "snv_fpr": r(sv["fpr"]),
        "indel_tp": ind["tp"], "indel_fp": ind["fp"], "indel_fn": ind["fn"],
        "indel_sensitivity": r(ind["sensitivity"]), "indel_precision": r(ind["precision"]),
        "indel_specificity": r(ind["specificity"]), "indel_fpr": r(ind["fpr"]),
        # ML-filtered metrics (NA when --ml-model not used)
        "ml_tp":  ml_ov["tp"]  if has_ml else "NA",
        "ml_fp":  ml_ov["fp"]  if has_ml else "NA",
        "ml_fn":  ml_ov["fn"]  if has_ml else "NA",
        "ml_sensitivity": r(ml_ov["sensitivity"]) if has_ml else "NA",
        "ml_precision":   r(ml_ov["precision"])   if has_ml else "NA",
        "ml_f1":          r(ml_ov["f1"])           if has_ml else "NA",
        "ml_snv_tp":  ml_sv["tp"]  if has_ml else "NA",
        "ml_snv_fp":  ml_sv["fp"]  if has_ml else "NA",
        "ml_snv_fn":  ml_sv["fn"]  if has_ml else "NA",
        "ml_snv_sensitivity": r(ml_sv["sensitivity"]) if has_ml else "NA",
        "ml_indel_tp":  ml_ind["tp"]  if has_ml else "NA",
        "ml_indel_fp":  ml_ind["fp"]  if has_ml else "NA",
        "ml_indel_fn":  ml_ind["fn"]  if has_ml else "NA",
        "ml_indel_sensitivity": r(ml_ind["sensitivity"]) if has_ml else "NA",
        "t_assemble_ms": stage_ms["assemble"],
        "t_index_ms": stage_ms["index"],
        "t_pathfind_ms": stage_ms["pathfind"],
        "t_call_ms": stage_ms["call"],
        "t_output_ms": stage_ms["output"],
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
        "peak_rss_mb": round(peak_mb, 1),
        "wall_time_s": round(elapsed, 1),
        "exit_code": proc.returncode,
    }

    if VCF_SAVE_DIR is not None and proc.returncode == 0:
        VCF_SAVE_DIR.mkdir(parents=True, exist_ok=True)
        src_vcf = out_dir / "variants.vcf"
        if src_vcf.exists():
            if KAM_TARGET_VARIANTS is not None:
                dst = VCF_SAVE_DIR / f"{name}.monitoring.vcf"
                shutil.copy(src_vcf, dst)
                disc_out = tmp_dir / "kam_disc"
                disc_out.mkdir()
                disc_cmd = [
                    str(KAM), "run",
                    "--r1", str(r1_sub), "--r2", str(r2_sub),
                    "--targets", str(TARGETS),
                    "--output-dir", str(disc_out),
                    "--output-format-override", "vcf",
                ]
                if KAM_MAX_VAF is not None:
                    disc_cmd += ["--max-vaf", str(KAM_MAX_VAF)]
                if KAM_MIN_ALT_MOLECULES is not None:
                    disc_cmd += ["--min-alt-molecules", str(KAM_MIN_ALT_MOLECULES)]
                if KAM_MIN_CONFIDENCE is not None:
                    disc_cmd += ["--min-confidence", str(KAM_MIN_CONFIDENCE)]
                if KAM_MIN_FAMILY_SIZE is not None:
                    disc_cmd += ["--min-family-size", str(KAM_MIN_FAMILY_SIZE)]
                if KAM_KMER_SIZE is not None:
                    disc_cmd += ["-k", str(KAM_KMER_SIZE)]
                if KAM_MIN_ALT_DUPLEX is not None:
                    disc_cmd += ["--min-alt-duplex", str(KAM_MIN_ALT_DUPLEX)]
                disc_proc = subprocess.run(disc_cmd, capture_output=True)
                disc_vcf = disc_out / "variants.vcf"
                if disc_vcf.exists():
                    shutil.copy(disc_vcf, VCF_SAVE_DIR / f"{name}.discovery.vcf")
            else:
                shutil.copy(src_vcf, VCF_SAVE_DIR / f"{name}.discovery.vcf")

    # ── Per-target TSV ───────────────────────────────────────────────────────
    if PER_SAMPLE_DIR is not None and proc.returncode == 0:
        per_target_path = PER_SAMPLE_DIR / f"{name}.targets.tsv"
        write_per_target_tsv(sample, scores, *truth_set, per_target_path)

    ml_sens_str = f" ml_sens={ml_ov['sensitivity']:.3f}" if has_ml else ""
    status = "OK" if proc.returncode == 0 else "FAIL"
    print(f"  [{name}] {status} | "
          f"mols={molecules:,} | "
          f"sens={ov['sensitivity']:.3f} spec={ov['specificity']:.3f} "
          f"prec={ov['precision']:.3f}{ml_sens_str} | "
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
    global KAM_MAX_VAF, KAM_MIN_ALT_MOLECULES, KAM_MIN_CONFIDENCE, KAM_MIN_FAMILY_SIZE
    global KAM_TARGET_VARIANTS, KAM_KMER_SIZE, KAM_MIN_ALT_DUPLEX, VCF_SAVE_DIR
    global KAM_ML_MODEL, PER_SAMPLE_DIR, KAM_TI_RESCUE

    parser = argparse.ArgumentParser(
        description="Run kam on all titration samples and score against truth variants."
    )
    parser.add_argument("--fastq-dir", type=Path, default=_DEFAULT_FASTQ)
    parser.add_argument("--truth-vcf", type=Path, default=_DEFAULT_TRUTH)
    parser.add_argument("--targets", type=Path, default=_DEFAULT_TARGETS)
    parser.add_argument("--kam-binary", type=Path, default=_DEFAULT_KAM)
    parser.add_argument("--results-dir", type=Path, default=_DEFAULT_RESULTS)
    parser.add_argument("--output", type=str, default=None)
    parser.add_argument("--reads", type=int, default=READS_PER_SAMPLE)
    parser.add_argument("--rss-limit-gb", type=float, default=PEAK_RSS_LIMIT_MB / 1024)
    parser.add_argument("--max-vaf", type=float, default=None)
    parser.add_argument("--min-alt-molecules", type=int, default=None)
    parser.add_argument("--min-confidence", type=float, default=None)
    parser.add_argument("--min-family-size", type=int, default=None)
    parser.add_argument("--target-variants", type=Path, default=None)
    parser.add_argument("--min-alt-duplex", type=int, default=None)
    parser.add_argument("--save-vcfs", type=Path, default=None)
    parser.add_argument("--kmer-size", type=int, default=None)
    parser.add_argument("--ml-model", type=str, default=None,
                        help="Built-in ML model name to pass to kam (e.g. twist-duplex-v1). "
                             "Adds ml_prob/ml_filter columns to output and computes ML-filtered metrics.")
    parser.add_argument("--per-sample-dir", type=Path, default=None,
                        help="Directory for per-sample per-target TSV files. "
                             "Each sample produces <name>.targets.tsv with one row per truth variant.")
    parser.add_argument("--ti-rescue", action="store_true", default=False,
                        help="Pass --ti-rescue to kam. Adds RESCUED/NO_EVIDENCE rows for undetected TI targets.")
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
    KAM_MIN_CONFIDENCE    = args.min_confidence
    KAM_MIN_FAMILY_SIZE   = args.min_family_size
    KAM_TARGET_VARIANTS   = args.target_variants
    KAM_MIN_ALT_DUPLEX    = args.min_alt_duplex
    KAM_KMER_SIZE         = args.kmer_size
    VCF_SAVE_DIR          = args.save_vcfs
    KAM_ML_MODEL          = args.ml_model
    PER_SAMPLE_DIR        = args.per_sample_dir
    KAM_TI_RESCUE         = args.ti_rescue

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
        "tp", "fp", "fn", "tn",
        "sensitivity", "precision", "specificity", "fpr", "fdr", "fnr", "f1",
        "snv_tp", "snv_fp", "snv_fn",
        "snv_sensitivity", "snv_precision", "snv_specificity", "snv_fpr",
        "indel_tp", "indel_fp", "indel_fn",
        "indel_sensitivity", "indel_precision", "indel_specificity", "indel_fpr",
        # ML-filtered metrics
        "ml_tp", "ml_fp", "ml_fn",
        "ml_sensitivity", "ml_precision", "ml_f1",
        "ml_snv_tp", "ml_snv_fp", "ml_snv_fn", "ml_snv_sensitivity",
        "ml_indel_tp", "ml_indel_fp", "ml_indel_fn", "ml_indel_sensitivity",
        "t_assemble_ms", "t_index_ms", "t_pathfind_ms", "t_call_ms", "t_output_ms",
        "rss_assemble_mb", "rss_index_mb", "rss_pathfind_mb", "rss_call_mb", "rss_output_mb",
        "cpu_assemble_pct", "cpu_index_pct", "cpu_pathfind_pct", "cpu_call_pct", "cpu_output_pct",
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
