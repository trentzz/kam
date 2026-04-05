"""Run the full pipeline for the ML3 test split.

Steps:
  1. Simulate: run varforge on all 1000 test configs.
  2. Kam: run kam discovery + tumour-informed on each simulation.
  3. Samples: build per-sample directories with truth TSV, calls TSVs, params.json.

Output directories:
  bigdata/experiments/02-ml-single-strand/results/test/sim_<name>/   — varforge outputs
  bigdata/experiments/02-ml-single-strand/results/test/kam_<name>/   — kam outputs
  bigdata/experiments/02-ml-single-strand/samples/test/<name>/        — sample dirs for training data builder

Usage: python3 scripts/ml/run_ml3_test_pipeline.py
Run from the repository root.
"""

import gzip
import json
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

REPO    = Path("/home/trent/code/kam")
CONFIGS = sorted((REPO / "bigdata/experiments/02-ml-single-strand/configs/test").glob("*.yaml"))
SIM_DIR = REPO / "bigdata/experiments/02-ml-single-strand/results/test"
KAM_DIR = REPO / "bigdata/experiments/02-ml-single-strand/results/test"
SAMPLES = REPO / "bigdata/experiments/02-ml-single-strand/samples/test"
KAM     = REPO / "target/release/kam"
VARFORGE = shutil.which("varforge") or "varforge"

TARGETS = {
    "snv":    REPO / "docs/benchmarking/01-snvindel/data/snvindel_targets.fa",
    "indel":  REPO / "docs/benchmarking/01-snvindel/data/snvindel_targets.fa",
    "sv":     REPO / "docs/benchmarking/02-sv-core/data/sv_suite_targets.fa",
    "ins":    REPO / "docs/benchmarking/02-sv-core/data/ins_targets.fa",
    "invdel": REPO / "docs/benchmarking/02-sv-core/data/invdel_targets.fa",
}


def variant_type(name: str) -> str:
    for vt in ("snv", "indel", "sv", "ins", "invdel"):
        if name.startswith(vt + "_"):
            return vt
    return "snv"


# ── Step 1: varforge simulation ────────────────────────────────────────────────

def simulate_one(config: Path) -> tuple[str, str]:
    name = config.stem
    out_dir = SIM_DIR / f"sim_{name}"
    r1 = out_dir / f"{name.upper()}_R1.fastq.gz"
    if r1.exists():
        return name, "skip"
    out_dir.mkdir(parents=True, exist_ok=True)
    result = subprocess.run(
        [VARFORGE, "simulate", "--config", str(config)],
        capture_output=True, text=True, cwd=REPO,
    )
    if result.returncode != 0:
        return name, f"sim_fail:{result.stderr[-300:]}"
    return name, "ok"


def step1_simulate(workers: int = 6) -> int:
    print("\n=== Step 1: varforge simulation ===", flush=True)
    to_run = [c for c in CONFIGS if not (SIM_DIR / f"sim_{c.stem}" / f"{c.stem.upper()}_R1.fastq.gz").exists()]
    print(f"Running {len(to_run)} simulations ({len(CONFIGS) - len(to_run)} already done)...", flush=True)
    done = skipped = failed = 0
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(simulate_one, c): c for c in to_run}
        for fut in as_completed(futures):
            name, status = fut.result()
            if status == "skip":
                skipped += 1
            elif status == "ok":
                done += 1
                print(f"[SIM] {name}", flush=True)
            else:
                failed += 1
                print(f"[SIM FAIL] {name}: {status}", flush=True)
    print(f"Step 1 done. Simulated={done} Skipped={skipped} Failed={failed}", flush=True)
    return failed


# ── Step 2: kam pipeline ───────────────────────────────────────────────────────

def run_kam_one(name: str) -> tuple[str, str]:
    vtype = variant_type(name)
    targets = TARGETS.get(vtype, TARGETS["snv"])
    sim_out = SIM_DIR / f"sim_{name}"
    r1 = sim_out / f"{name.upper()}_R1.fastq.gz"
    r2 = sim_out / f"{name.upper()}_R2.fastq.gz"

    if not r1.exists():
        return name, "no_fastq"

    kam_out = KAM_DIR / f"kam_{name}"
    kam_out.mkdir(parents=True, exist_ok=True)

    # Discovery run.
    disc_dir = kam_out / "tmp_disc"
    disc_tsv = kam_out / "calls_discovery.tsv"
    disc_log = kam_out / "discovery.log"
    if not disc_tsv.exists():
        disc_dir.mkdir(parents=True, exist_ok=True)
        r = subprocess.run(
            [str(KAM), "run",
             "--r1", str(r1), "--r2", str(r2),
             "--targets", str(targets),
             "--output-dir", str(disc_dir),
             "--output-format-override", "vcf,tsv"],
            capture_output=True, text=True, cwd=REPO,
        )
        # Save full stderr log (QC progress lines from kam).
        disc_log.write_text(r.stderr)
        if r.returncode != 0:
            return name, f"disc_fail:{r.stderr[-400:]}"
        # Copy all output files (TSV, VCF, QC JSONs) to kam_out.
        for f in disc_dir.iterdir():
            shutil.copy(f, kam_out / f"discovery_{f.name}")
        for f in disc_dir.glob("*.tsv"):
            shutil.copy(f, disc_tsv)

    # Decompress truth VCF.
    truth_vcf_gz = sim_out / f"{name.upper()}.truth.vcf.gz"
    truth_vcf    = sim_out / f"{name.upper()}.truth.vcf"
    if truth_vcf_gz.exists() and not truth_vcf.exists():
        with gzip.open(truth_vcf_gz, "rb") as fin, open(truth_vcf, "wb") as fout:
            shutil.copyfileobj(fin, fout)

    # Tumour-informed run.
    ti_dir = kam_out / "tmp_ti"
    ti_tsv = kam_out / "calls_tumour_informed.tsv"
    ti_log = kam_out / "tumour_informed.log"
    if truth_vcf.exists() and not ti_tsv.exists():
        ti_dir.mkdir(parents=True, exist_ok=True)
        r = subprocess.run(
            [str(KAM), "run",
             "--r1", str(r1), "--r2", str(r2),
             "--targets", str(targets),
             "--output-dir", str(ti_dir),
             "--output-format-override", "vcf,tsv",
             "--target-variants", str(truth_vcf)],
            capture_output=True, text=True, cwd=REPO,
        )
        # Save full stderr log.
        ti_log.write_text(r.stderr)
        if r.returncode != 0:
            return name, f"ti_fail:{r.stderr[-400:]}"
        # Copy all output files (TSV, VCF, QC JSONs) to kam_out.
        for f in ti_dir.iterdir():
            shutil.copy(f, kam_out / f"tumour_informed_{f.name}")
        for f in ti_dir.glob("*.tsv"):
            shutil.copy(f, ti_tsv)

    return name, "ok"


def step2_kam(workers: int = 12) -> int:
    print("\n=== Step 2: kam runs ===", flush=True)
    names_with_fastq = [
        c.stem for c in CONFIGS
        if (SIM_DIR / f"sim_{c.stem}" / f"{c.stem.upper()}_R1.fastq.gz").exists()
    ]
    to_run = [n for n in names_with_fastq if not (KAM_DIR / f"kam_{n}" / "calls_discovery.tsv").exists()]
    print(f"Running kam on {len(to_run)} samples ({len(names_with_fastq) - len(to_run)} already done)...", flush=True)
    done = skipped = failed = 0
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(run_kam_one, n): n for n in to_run}
        for fut in as_completed(futures):
            name, status = fut.result()
            if status == "ok":
                done += 1
                print(f"[KAM] {name}", flush=True)
            elif status == "no_fastq":
                skipped += 1
            else:
                failed += 1
                print(f"[KAM FAIL] {name}: {status}", flush=True)
    print(f"Step 2 done. Ran={done} Skipped={skipped} Failed={failed}", flush=True)
    return failed


# ── Step 3: build sample directories ──────────────────────────────────────────

def build_sample_one(name: str) -> tuple[str, str]:
    vtype  = variant_type(name)
    sim_out = SIM_DIR / f"sim_{name}"
    kam_out = KAM_DIR / f"kam_{name}"
    sample_dir = SAMPLES / name

    disc_tsv = kam_out / "calls_discovery.tsv"
    if not disc_tsv.exists():
        return name, "no_calls"

    sample_dir.mkdir(parents=True, exist_ok=True)

    # Truth TSV from truth VCF.
    truth_vcf = sim_out / f"{name.upper()}.truth.vcf"
    if not truth_vcf.exists():
        truth_vcf_gz = sim_out / f"{name.upper()}.truth.vcf.gz"
        if truth_vcf_gz.exists():
            with gzip.open(truth_vcf_gz, "rb") as fin, open(truth_vcf, "wb") as fout:
                shutil.copyfileobj(fin, fout)

    if truth_vcf.exists():
        truth_tsv = sample_dir / "truth.tsv"
        with open(truth_vcf) as f:
            lines = [l for l in f if not l.startswith("#")]
        with open(truth_tsv, "w") as f:
            f.write("chrom\tpos\tref\talt\n")
            for line in lines:
                parts = line.strip().split("\t")
                if len(parts) >= 5:
                    chrom, pos, _, ref, alt = parts[:5]
                    f.write(f"{chrom}\t{int(pos) - 1}\t{ref}\t{alt}\n")

    # Copy TSVs.
    for src, dst in [
        (disc_tsv, sample_dir / "discovery.tsv"),
        (kam_out / "calls_tumour_informed.tsv", sample_dir / "tumour_informed.tsv"),
    ]:
        if src.exists():
            shutil.copy(src, dst)

    # Rename columns in TSVs to match build_training_data_v2.py expectations.
    for tsv in [sample_dir / "discovery.tsv", sample_dir / "tumour_informed.tsv"]:
        if not tsv.exists():
            continue
        text = tsv.read_text()
        text = (text
                .replace("target_id", "chrom")
                .replace("vaf_ci_low", "vaf_lo")
                .replace("vaf_ci_high", "vaf_hi")
                .replace("n_molecules_ref", "nref")
                .replace("n_molecules_alt", "nalt")
                .replace("n_duplex_alt", "ndupalt")
                .replace("n_simplex_alt", "nsimalt")
                .replace("strand_bias_p", "sbp")
                .replace("confidence", "conf")
                .replace("ref_seq", "ref")
                .replace("alt_seq", "alt"))
        tsv.write_text(text)

    # Parse position from chrom column (chr:start-end → start).
    for tsv in [sample_dir / "discovery.tsv", sample_dir / "tumour_informed.tsv"]:
        if not tsv.exists():
            continue
        import csv, io
        rows = list(csv.DictReader(tsv.read_text(), delimiter="\t"))
        if not rows or "chrom" not in rows[0]:
            continue
        for row in rows:
            chrom = row.get("chrom", "")
            if ":" in chrom:
                parts = chrom.split(":")
                row["chrom"] = parts[0]
                coords = parts[1].split("-")
                row["pos"] = coords[0] if coords else "0"
        out = io.StringIO()
        w = csv.DictWriter(out, fieldnames=rows[0].keys(), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
        tsv.write_text(out.getvalue())

    # Copy params.json.
    params_src = REPO / "bigdata/experiments/02-ml-single-strand/configs/test" / f"{name}_params.json"
    if params_src.exists():
        shutil.copy(params_src, sample_dir / "params.json")

    return name, "ok"


def step3_samples() -> int:
    print("\n=== Step 3: build sample directories ===", flush=True)
    names = [c.stem for c in CONFIGS if (KAM_DIR / f"kam_{c.stem}" / "calls_discovery.tsv").exists()]
    to_build = [n for n in names if not (SAMPLES / n / "discovery.tsv").exists()]
    print(f"Building {len(to_build)} sample dirs ({len(names) - len(to_build)} already done)...", flush=True)
    done = skipped = failed = 0
    for name in to_build:
        n, status = build_sample_one(name)
        if status == "ok":
            done += 1
            print(f"[SAMPLE] {n}", flush=True)
        else:
            failed += 1
            print(f"[SAMPLE FAIL] {n}: {status}", flush=True)
    print(f"Step 3 done. Built={done} Skipped={skipped} Failed={failed}", flush=True)
    return failed


if __name__ == "__main__":
    SIM_DIR.mkdir(parents=True, exist_ok=True)
    SAMPLES.mkdir(parents=True, exist_ok=True)

    step1_simulate(workers=6)
    step2_kam(workers=12)
    step3_samples()

    print("\n=== ML3 test pipeline complete ===", flush=True)
