"""Full ML data generation pipeline for kam.

Steps:
  1. Run varforge on all configs that don't have FASTQs yet.
  2. Run kam (discovery + tumour-informed) on all sims with FASTQs.
  3. Build per-sample directories with truth TSV, calls TSVs, and params.json.
  4. Run build_training_data_v2.py to generate training_data_v2.csv.
"""

import gzip
import json
import re
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

REPO = Path("/home/trent/code/kam")
CONFIGS = sorted((REPO / "docs/benchmarking/ml/configs").glob("*.yaml"))
RESULTS = REPO / "docs/benchmarking/ml/results"
SAMPLES = REPO / "docs/benchmarking/ml/samples"
KAM = REPO / "target/release/kam"

TARGETS = {
    "snv": REPO / "docs/benchmarking/snvindel/data/snvindel_targets.fa",
    "indel": REPO / "docs/benchmarking/snvindel/data/snvindel_targets.fa",
    "sv": REPO / "docs/benchmarking/sv/data/sv_suite_targets.fa",
    "ins": REPO / "docs/benchmarking/sv/data/ins_targets.fa",
    "invdel": REPO / "docs/benchmarking/sv/data/invdel_targets.fa",
}


# ---------------------------------------------------------------------------
# Step 1: varforge simulation
# ---------------------------------------------------------------------------

def needs_sim(config_path):
    name = config_path.stem
    sim_dir = RESULTS / f"sim_{name}"
    return not any(sim_dir.glob("*_R1.fastq.gz"))


def run_sim(config_path):
    name = config_path.stem
    sim_dir = RESULTS / f"sim_{name}"
    if not needs_sim(config_path):
        return name, "skipped"
    sim_dir.mkdir(parents=True, exist_ok=True)
    result = subprocess.run(
        ["varforge", "simulate", "--config", str(config_path.relative_to(REPO))],
        cwd=str(REPO), capture_output=True, text=True
    )
    if result.returncode != 0:
        print(f"[SIM FAIL] {name}: {result.stderr[-400:]}", flush=True)
        return name, "failed"
    print(f"[SIM] {name}", flush=True)
    return name, "done"


def step1_simulate(workers=6):
    print("\n=== Step 1: varforge simulation ===", flush=True)
    to_run = [c for c in CONFIGS if needs_sim(c)]
    print(f"Running {len(to_run)} simulations ({len(CONFIGS) - len(to_run)} already done)...", flush=True)
    done = skipped = failed = 0
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(run_sim, c): c for c in CONFIGS}
        for fut in as_completed(futures):
            name, status = fut.result()
            if status == "done":
                done += 1
            elif status == "skipped":
                skipped += 1
            else:
                failed += 1
    print(f"Step 1 done. Simulated={done} Skipped={skipped} Failed={failed}", flush=True)
    return failed


# ---------------------------------------------------------------------------
# Step 2: kam runs
# ---------------------------------------------------------------------------

def get_variant_type(name):
    """Parse variant type from sample name prefix."""
    for t in ("snv", "indel", "sv", "ins", "invdel"):
        if name.startswith(t + "_"):
            return t
    # Try after removing leading underscores or dashes
    for t in ("invdel", "ins", "sv", "indel", "snv"):
        if t in name:
            return t
    return None


def find_fastqs(sim_dir):
    r1s = sorted(sim_dir.glob("*_R1.fastq.gz"))
    r2s = sorted(sim_dir.glob("*_R2.fastq.gz"))
    if not r1s or not r2s:
        return None, None
    return r1s[0], r2s[0]


def find_truth_vcf(sim_dir):
    """Return path to truth VCF, preferring uncompressed.

    If only a .gz exists, decompress it alongside the .gz file and return
    the plain path. kam's --target-variants flag requires an uncompressed VCF.
    """
    vcf = sorted(sim_dir.glob("*.truth.vcf"))
    if vcf:
        return vcf[0]
    vcfgz = sorted(sim_dir.glob("*.truth.vcf.gz"))
    if vcfgz:
        # Decompress to sibling file
        plain = vcfgz[0].with_suffix("")  # strips .gz
        if not plain.exists():
            with gzip.open(vcfgz[0], "rt") as fin, open(plain, "w") as fout:
                fout.write(fin.read())
        return plain
    return None


def needs_kam(name):
    kam_dir = RESULTS / f"kam_{name}"
    return not (kam_dir / "calls_discovery.tsv").exists()


def run_kam(name):
    if not needs_kam(name):
        return name, "skipped"

    sim_dir = RESULTS / f"sim_{name}"
    r1, r2 = find_fastqs(sim_dir)
    if not r1:
        return name, "no_fastq"

    truth_vcf = find_truth_vcf(sim_dir)
    if not truth_vcf:
        return name, "no_truth_vcf"

    vtype = get_variant_type(name)
    if vtype not in TARGETS:
        print(f"[KAM WARN] Unknown variant type for {name}", flush=True)
        return name, "unknown_type"

    targets = TARGETS[vtype]
    kam_dir = RESULTS / f"kam_{name}"
    kam_dir.mkdir(parents=True, exist_ok=True)

    disc_tmp = kam_dir / "tmp_disc"
    ti_tmp = kam_dir / "tmp_ti"

    # Discovery run
    disc_result = subprocess.run(
        [str(KAM), "run",
         "--r1", str(r1), "--r2", str(r2),
         "--targets", str(targets),
         "--output-dir", str(disc_tmp),
         "--output-format-override", "vcf,tsv"],
        cwd=str(REPO), capture_output=True, text=True
    )
    if disc_result.returncode != 0:
        print(f"[KAM DISC FAIL] {name}: {disc_result.stderr[-400:]}", flush=True)
        return name, "disc_failed"

    # Tumour-informed run
    ti_result = subprocess.run(
        [str(KAM), "run",
         "--r1", str(r1), "--r2", str(r2),
         "--targets", str(targets),
         "--target-variants", str(truth_vcf),
         "--output-dir", str(ti_tmp),
         "--output-format-override", "vcf,tsv"],
        cwd=str(REPO), capture_output=True, text=True
    )
    if ti_result.returncode != 0:
        print(f"[KAM TI FAIL] {name}: {ti_result.stderr[-400:]}", flush=True)
        return name, "ti_failed"

    # Copy outputs
    for src, dst in [
        (disc_tmp / "variants.tsv", kam_dir / "calls_discovery.tsv"),
        (disc_tmp / "variants.vcf", kam_dir / "calls_discovery.vcf"),
        (ti_tmp / "variants.tsv", kam_dir / "calls_tumour_informed.tsv"),
        (ti_tmp / "variants.vcf", kam_dir / "calls_tumour_informed.vcf"),
    ]:
        if src.exists():
            shutil.copy2(src, dst)

    # Remove tmp dirs
    shutil.rmtree(disc_tmp, ignore_errors=True)
    shutil.rmtree(ti_tmp, ignore_errors=True)

    print(f"[KAM] {name}", flush=True)
    return name, "done"


def step2_kam(workers=4):
    print("\n=== Step 2: kam runs ===", flush=True)
    # Only run on sims that have FASTQs
    sim_dirs = [d for d in RESULTS.iterdir() if d.is_dir() and d.name.startswith("sim_")]
    names_with_fastq = [d.name[4:] for d in sim_dirs if any(d.glob("*_R1.fastq.gz"))]
    to_run = [n for n in names_with_fastq if needs_kam(n)]
    print(f"Running kam on {len(to_run)} samples ({len(names_with_fastq) - len(to_run)} already done)...", flush=True)
    done = skipped = failed = 0
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(run_kam, n): n for n in names_with_fastq}
        for fut in as_completed(futures):
            name, status = fut.result()
            if status == "done":
                done += 1
            elif status == "skipped":
                skipped += 1
            else:
                failed += 1
    print(f"Step 2 done. Ran={done} Skipped={skipped} Failed={failed}", flush=True)
    return failed


# ---------------------------------------------------------------------------
# Step 3: build sample directories
# ---------------------------------------------------------------------------

def parse_truth_vcf(vcf_path):
    """Parse varforge truth VCF to list of dicts."""
    rows = []
    opener = gzip.open if str(vcf_path).endswith(".gz") else open
    with opener(vcf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip().split("\t")
            chrom, pos, _, ref, alt = parts[0], parts[1], parts[2], parts[3], parts[4]
            info = parts[7]
            # Parse EXPECTED_VAF or VAF
            vaf_match = re.search(r"EXPECTED_VAF=([0-9.eE+\-]+)", info)
            if not vaf_match:
                vaf_match = re.search(r"VAF=([0-9.eE+\-]+)", info)
            vaf = float(vaf_match.group(1)) if vaf_match else ""
            # Parse VARTYPE or TYPE
            type_match = re.search(r"VARTYPE=([^;]+)", info)
            if not type_match:
                type_match = re.search(r"TYPE=([^;]+)", info)
            vartype = type_match.group(1) if type_match else ""
            rows.append({"chrom": chrom, "pos": pos, "ref": ref, "alt": alt, "vaf": vaf, "type": vartype})
    return rows


def parse_config_params(config_path, name):
    """Parse params.json fields from config YAML (without yaml dep)."""
    content = config_path.read_text()
    params = {}

    # coverage
    cov_m = re.search(r"^\s*coverage:\s*([\d.]+)", content, re.MULTILINE)
    params["coverage"] = float(cov_m.group(1)) if cov_m else None

    # family_size_mean
    fsm = re.search(r"^\s*family_size_mean:\s*([\d.]+)", content, re.MULTILINE)
    params["family_size_mean"] = float(fsm.group(1)) if fsm else None

    # pcr_cycles
    pcr = re.search(r"^\s*pcr_cycles:\s*([\d.]+)", content, re.MULTILINE)
    params["pcr_cycles"] = int(float(pcr.group(1))) if pcr else None

    # fragment mean
    frag = re.search(r"^\s*mean:\s*([\d.]+)", content, re.MULTILINE)
    params["fragment_mean"] = float(frag.group(1)) if frag else None

    # Parse variant type and VAF from name
    # name format: snv_vaf0100_c or snv_vaf0100_cov1000_a etc.
    # vtype is the prefix
    vtype = get_variant_type(name)
    params["variant_type"] = vtype

    # VAF: from vafXXXX in name (integer, divide by 10000 for fraction)
    vaf_m = re.search(r"vaf(\d+)", name)
    if vaf_m:
        vaf_int = int(vaf_m.group(1))
        params["vaf"] = vaf_int / 10000.0
    else:
        params["vaf"] = None

    # Replicate: last character of the name (a, b, c, d, e, ...)
    # But param_variant may change the "base" key
    # Detect param variant from name parts after vafXXXX
    # e.g. snv_vaf0100_c -> replicate=c, param_variant=base
    # snv_vaf0100_cov1000_a -> replicate=a, param_variant=cov1000
    parts_after_vaf = name.split(f"vaf{vaf_m.group(1)}_")[1] if vaf_m else ""
    # param_variant keywords
    param_variant_match = re.match(r"(cov\d+|fam\d+|pcr\d+|frag\d+)_(.+)", parts_after_vaf)
    if param_variant_match:
        params["param_variant"] = param_variant_match.group(1)
        params["replicate"] = param_variant_match.group(2)
    else:
        # Just a replicate letter
        params["param_variant"] = "base"
        params["replicate"] = parts_after_vaf if parts_after_vaf else "a"

    return params


def build_sample_dir(name, config_path):
    sim_dir = RESULTS / f"sim_{name}"
    kam_dir = RESULTS / f"kam_{name}"
    sample_dir = SAMPLES / name

    # Check prerequisites
    if not (kam_dir / "calls_discovery.tsv").exists():
        return name, "no_kam"

    truth_vcf = find_truth_vcf(sim_dir)
    if not truth_vcf:
        return name, "no_truth_vcf"

    sample_dir.mkdir(parents=True, exist_ok=True)

    # config.yaml
    shutil.copy2(config_path, sample_dir / "config.yaml")

    # varforge_cmd.txt
    cmd_path = config_path.parent / f"{name}_cmd.txt"
    if cmd_path.exists():
        shutil.copy2(cmd_path, sample_dir / "varforge_cmd.txt")
    else:
        (sample_dir / "varforge_cmd.txt").write_text(
            f"varforge simulate --config docs/benchmarking/ml/configs/{name}.yaml\n"
        )

    # truth.tsv
    rows = parse_truth_vcf(truth_vcf)
    with open(sample_dir / "truth.tsv", "w") as fh:
        fh.write("chrom\tpos\tref\talt\tvaf\ttype\n")
        for r in rows:
            fh.write(f"{r['chrom']}\t{r['pos']}\t{r['ref']}\t{r['alt']}\t{r['vaf']}\t{r['type']}\n")

    # discovery.tsv and tumour_informed.tsv
    shutil.copy2(kam_dir / "calls_discovery.tsv", sample_dir / "discovery.tsv")
    ti_src = kam_dir / "calls_tumour_informed.tsv"
    if ti_src.exists():
        shutil.copy2(ti_src, sample_dir / "tumour_informed.tsv")

    # params.json
    params = parse_config_params(config_path, name)
    with open(sample_dir / "params.json", "w") as fh:
        json.dump(params, fh, indent=2)

    print(f"[SAMPLE] {name}", flush=True)
    return name, "done"


def step3_samples():
    print("\n=== Step 3: build sample directories ===", flush=True)
    SAMPLES.mkdir(parents=True, exist_ok=True)

    # Map config stem to path
    config_map = {c.stem: c for c in CONFIGS}

    done = skipped = failed = 0
    for name, config_path in config_map.items():
        kam_dir = RESULTS / f"kam_{name}"
        if not (kam_dir / "calls_discovery.tsv").exists():
            skipped += 1
            continue
        _, status = build_sample_dir(name, config_path)
        if status == "done":
            done += 1
        else:
            print(f"[SAMPLE FAIL] {name}: {status}", flush=True)
            failed += 1

    print(f"Step 3 done. Built={done} Skipped={skipped} Failed={failed}", flush=True)
    return failed


# ---------------------------------------------------------------------------
# Step 4: build_training_data_v2.py
# ---------------------------------------------------------------------------

def step4_training_data():
    print("\n=== Step 4: build_training_data_v2.py ===", flush=True)
    result = subprocess.run(
        ["python3", "scripts/ml/build_training_data_v2.py"],
        cwd=str(REPO), capture_output=False, text=True
    )
    if result.returncode != 0:
        print(f"[TRAINING DATA FAIL] returncode={result.returncode}", flush=True)
        return 1
    print("Step 4 done.", flush=True)
    return 0


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import sys

    sim_failures = step1_simulate(workers=6)
    kam_failures = step2_kam(workers=4)
    sample_failures = step3_samples()
    training_failure = step4_training_data()

    print(f"\n=== Pipeline complete ===", flush=True)
    print(f"Sim failures: {sim_failures}", flush=True)
    print(f"Kam failures: {kam_failures}", flush=True)
    print(f"Sample failures: {sample_failures}", flush=True)
    print(f"Training data failure: {training_failure}", flush=True)

    if any([sim_failures, kam_failures, sample_failures, training_failure]):
        sys.exit(1)
