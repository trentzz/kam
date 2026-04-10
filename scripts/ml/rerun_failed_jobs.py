"""Targeted re-run of 11 failed ML3 train pipeline jobs.

Category 1 (DFS-hang): delete stale tmp_disc/tmp_ti, then re-run kam.
Category 2 (invdel I/O errors): delete sim + kam dirs, re-simulate, then run kam.

Run from the repository root.
"""

import gzip
import json
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

REPO     = Path("/home/trent/code/kam")
SIM_DIR  = REPO / "bigdata/experiments/02-ml-single-strand/results/train"
KAM_DIR  = REPO / "bigdata/experiments/02-ml-single-strand/results/train"
SAMPLES  = REPO / "bigdata/experiments/02-ml-single-strand/samples/train"
KAM      = REPO / "target/release/kam"
VARFORGE = shutil.which("varforge") or "varforge"

TARGETS = {
    "snv":    REPO / "docs/benchmarking/01-snvindel/data/snvindel_targets.fa",
    "indel":  REPO / "docs/benchmarking/01-snvindel/data/snvindel_targets.fa",
    "sv":     REPO / "docs/benchmarking/02-sv-core/data/sv_suite_targets.fa",
    "ins":    REPO / "docs/benchmarking/02-sv-core/data/ins_targets.fa",
    "invdel": REPO / "docs/benchmarking/02-sv-core/data/invdel_targets.fa",
}

# Category 1: DFS-hang jobs — delete stale tmp dirs, re-run kam.
DFS_HANG_JOBS = [
    "ins_ml3_train_vaf0600_0953",
    "ins_ml3_train_vaf0807_1015",
    "ins_ml3_train_vaf1218_1534",
    "ins_ml3_train_vaf1380_0112",
    "ins_ml3_train_vaf0559_1678",
    "sv_ml3_train_vaf0469_1290",
]

# Category 2: invdel I/O error jobs — re-simulate then re-run kam.
INVDEL_JOBS = [
    "invdel_ml3_train_vaf0572_1213",
    "invdel_ml3_train_vaf0587_0280",
    "invdel_ml3_train_vaf0599_1664",
    "invdel_ml3_train_vaf0605_1465",
    "invdel_ml3_train_vaf0618_1427",
]


def variant_type(name: str) -> str:
    for vt in ("snv", "indel", "sv", "ins", "invdel"):
        if name.startswith(vt + "_"):
            return vt
    return "snv"


def run_kam_one(name: str) -> tuple[str, str]:
    """Run disc + tumour-informed for one job. Mirrors run_kam_one from the main pipeline."""
    vtype   = variant_type(name)
    targets = TARGETS.get(vtype, TARGETS["snv"])
    sim_out = SIM_DIR / f"sim_{name}"
    upper   = name.upper()
    r1      = sim_out / f"{upper}_R1.fastq.gz"
    r2      = sim_out / f"{upper}_R2.fastq.gz"

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
        disc_log.write_text(r.stderr)
        if r.returncode != 0:
            return name, f"disc_fail:{r.stderr[-400:]}"
        for f in disc_dir.iterdir():
            shutil.copy(f, kam_out / f"discovery_{f.name}")
        for f in disc_dir.glob("*.tsv"):
            shutil.copy(f, disc_tsv)

    # Decompress truth VCF.
    truth_vcf_gz = sim_out / f"{upper}.truth.vcf.gz"
    truth_vcf    = sim_out / f"{upper}.truth.vcf"
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
        ti_log.write_text(r.stderr)
        if r.returncode != 0:
            return name, f"ti_fail:{r.stderr[-400:]}"
        for f in ti_dir.iterdir():
            shutil.copy(f, kam_out / f"tumour_informed_{f.name}")
        for f in ti_dir.glob("*.tsv"):
            shutil.copy(f, ti_tsv)

    return name, "ok"


def handle_dfs_hang(name: str) -> tuple[str, str]:
    """Delete stale tmp dirs for a DFS-hang job, then re-run."""
    kam_out  = KAM_DIR / f"kam_{name}"
    disc_tsv = kam_out / "calls_discovery.tsv"
    ti_tsv   = kam_out / "calls_tumour_informed.tsv"

    # Delete stale tmp_disc only if discovery hasn't succeeded yet.
    if not disc_tsv.exists():
        tmp_disc = kam_out / "tmp_disc"
        if tmp_disc.exists():
            shutil.rmtree(tmp_disc)
            print(f"[CLEAN] Deleted {tmp_disc}", flush=True)

    # Delete stale tmp_ti only if tumour-informed hasn't succeeded yet.
    if not ti_tsv.exists():
        tmp_ti = kam_out / "tmp_ti"
        if tmp_ti.exists():
            shutil.rmtree(tmp_ti)
            print(f"[CLEAN] Deleted {tmp_ti}", flush=True)

    return run_kam_one(name)


def handle_invdel(name: str) -> tuple[str, str]:
    """Delete corrupt sim + kam dirs, re-simulate, then run kam."""
    sim_out = SIM_DIR / f"sim_{name}"
    kam_out = KAM_DIR / f"kam_{name}"
    cfg     = REPO / "bigdata/experiments/02-ml-single-strand/configs/train" / f"{name}.yaml"

    # Delete corrupt outputs.
    if sim_out.exists():
        shutil.rmtree(sim_out)
        print(f"[CLEAN] Deleted {sim_out}", flush=True)
    if kam_out.exists():
        shutil.rmtree(kam_out)
        print(f"[CLEAN] Deleted {kam_out}", flush=True)

    # Re-simulate.
    sim_out.mkdir(parents=True, exist_ok=True)
    result = subprocess.run(
        [VARFORGE, "simulate", "--config", str(cfg)],
        capture_output=True, text=True, cwd=REPO,
    )
    if result.returncode != 0:
        return name, f"sim_fail:{result.stderr[-400:]}"
    print(f"[SIM] {name}", flush=True)

    return run_kam_one(name)


def build_sample_one(name: str) -> tuple[str, str]:
    """Build the per-sample directory. Mirrors build_sample_one from the main pipeline."""
    vtype      = variant_type(name)
    sim_out    = SIM_DIR / f"sim_{name}"
    kam_out    = KAM_DIR / f"kam_{name}"
    sample_dir = SAMPLES / name

    disc_tsv = kam_out / "calls_discovery.tsv"
    if not disc_tsv.exists():
        return name, "no_calls"

    sample_dir.mkdir(parents=True, exist_ok=True)

    # Truth TSV from truth VCF.
    upper    = name.upper()
    truth_vcf = sim_out / f"{upper}.truth.vcf"
    if not truth_vcf.exists():
        truth_vcf_gz = sim_out / f"{upper}.truth.vcf.gz"
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

    # Rename columns to match build_training_data_v3.py expectations.
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

    # Copy params.json.
    params_src = REPO / "bigdata/experiments/02-ml-single-strand/configs/train" / f"{name}_params.json"
    if params_src.exists():
        shutil.copy(params_src, sample_dir / "params.json")

    return name, "ok"


if __name__ == "__main__":
    all_jobs = DFS_HANG_JOBS + INVDEL_JOBS

    print(f"=== Targeted re-run: {len(all_jobs)} jobs ===", flush=True)
    print(f"  Category 1 (DFS-hang): {DFS_HANG_JOBS}", flush=True)
    print(f"  Category 2 (invdel):   {INVDEL_JOBS}", flush=True)

    results: dict[str, str] = {}

    def run_job(name: str) -> tuple[str, str]:
        if name in INVDEL_JOBS:
            return handle_invdel(name)
        else:
            return handle_dfs_hang(name)

    with ThreadPoolExecutor(max_workers=6) as pool:
        futures = {pool.submit(run_job, n): n for n in all_jobs}
        for fut in as_completed(futures):
            name, status = fut.result()
            results[name] = status
            if status == "ok":
                print(f"[KAM OK]  {name}", flush=True)
            else:
                print(f"[KAM FAIL] {name}: {status}", flush=True)

    print("\n=== Building sample directories ===", flush=True)
    sample_results: dict[str, str] = {}
    for name in all_jobs:
        if results.get(name) == "ok":
            _, status = build_sample_one(name)
            sample_results[name] = status
            if status == "ok":
                print(f"[SAMPLE OK]   {name}", flush=True)
            else:
                print(f"[SAMPLE FAIL] {name}: {status}", flush=True)

    print("\n=== Summary ===", flush=True)
    ok_kam     = [n for n, s in results.items() if s == "ok"]
    fail_kam   = {n: s for n, s in results.items() if s != "ok"}
    ok_sample  = [n for n, s in sample_results.items() if s == "ok"]
    fail_sample = {n: s for n, s in sample_results.items() if s != "ok"}
    print(f"KAM:    OK={len(ok_kam)}  FAIL={len(fail_kam)}", flush=True)
    print(f"SAMPLE: OK={len(ok_sample)}  FAIL={len(fail_sample)}", flush=True)
    if fail_kam:
        print(f"KAM failures: {fail_kam}", flush=True)
    if fail_sample:
        print(f"SAMPLE failures: {fail_sample}", flush=True)
