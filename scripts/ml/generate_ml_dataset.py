#!/usr/bin/env python3
"""Generate a diverse ML training dataset using varforge and kam.

Produces many simulated samples with varied coverage, UMI family size, PCR
cycles, fragment size, and replicates. All output goes to
docs/benchmarking/ml/.

Usage:
    python3 scripts/ml/generate_ml_dataset.py

Run from the repository root.
"""

from __future__ import annotations

import gzip
import json
import os
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any

REPO = Path(__file__).resolve().parent.parent.parent

# ─── Paths ────────────────────────────────────────────────────────────────────

ML_DIR      = REPO / "docs" / "benchmarking" / "ml"
CONFIGS_DIR = ML_DIR / "configs"
RESULTS_DIR = ML_DIR / "results"
SAMPLES_DIR = ML_DIR / "samples"

SNVINDEL_DATA = REPO / "docs" / "benchmarking" / "snvindel" / "data"
SV_DATA       = REPO / "docs" / "benchmarking" / "sv" / "data"

SNVINDEL_TARGETS = SNVINDEL_DATA / "snvindel_targets.fa"
SV_TARGETS       = SV_DATA / "sv_suite_targets.fa"
INS_TARGETS      = SV_DATA / "ins_targets.fa"
INVDEL_TARGETS   = SV_DATA / "invdel_targets.fa"

KAM = REPO / "target" / "release" / "kam"

# ─── Base seeds ───────────────────────────────────────────────────────────────
# Formula: base_seed + rep_offset * 1000 + vaf_tag_int
# rep a=0, b=1, c=2, d=3, e=4

BASE_SEEDS = {
    "snv":    10000,
    "indel":  30000,
    "sv":     50000,
    "ins":    70000,
    "invdel": 90000,
}

REP_OFFSETS = {"a": 0, "b": 1, "c": 2, "d": 3, "e": 4}

# ─── VAF levels ───────────────────────────────────────────────────────────────

ALL_VAFS = [
    0.0005, 0.0010, 0.0020, 0.0030, 0.0050, 0.0075, 0.0100, 0.0150,
    0.0200, 0.0300, 0.0500, 0.0700, 0.1000, 0.0400, 0.0600,
]

SWEEP_VAFS = [
    0.0005, 0.0010, 0.0020, 0.0050, 0.0100, 0.0200, 0.0500, 0.1000,
]

# ─── Default params ───────────────────────────────────────────────────────────

DEFAULTS = {
    "coverage":         5000.0,
    "family_size_mean": 4.0,
    "family_size_sd":   1.5,
    "pcr_cycles":       8,
    "fragment_mean":    167.0,
    "fragment_sd":      30.0,
}

# ─── Helpers ──────────────────────────────────────────────────────────────────

def vaf_tag(vaf: float) -> str:
    """Convert a VAF float to a 4-digit tag string."""
    return f"{round(vaf * 10000):04d}"


def seed_for(vtype: str, rep: str, vaf: float) -> int:
    """Compute a deterministic seed from type, replicate, and VAF."""
    return BASE_SEEDS[vtype] + REP_OFFSETS[rep] * 1000 + round(vaf * 10000)


def truth_vcf_snvindel(vtype: str, tag: str, rep: str) -> Path:
    """Return the truth VCF path for an SNV/indel sample.

    For new replicates (c, d, e) that have no separate truth VCF, reuse rep=a.
    The seed change is what varies the stochastic output; positions are identical.
    """
    prefix = "truth_snvs" if vtype == "snv" else "truth_indels"
    # Replicates c/d/e reuse rep=a truth positions.
    base_rep = rep if rep in ("a", "b") else "a"
    return SNVINDEL_DATA / f"{prefix}_vaf{tag}_{base_rep}.vcf"


def truth_vcf_sv(vtype: str, tag: str, rep: str) -> Path:
    """Return the truth VCF path for an SV sample.

    Replicates c/d/e reuse rep=a positions.
    """
    prefixes = {"sv": "truth_svs", "ins": "truth_ins", "invdel": "truth_invdel"}
    base_rep = rep if rep in ("a", "b") else "a"
    return SV_DATA / f"{prefixes[vtype]}_vaf{tag}_{base_rep}.vcf"


def targets_for(vtype: str) -> Path:
    """Return the targets FASTA for a variant type."""
    if vtype in ("sv",):
        return SV_TARGETS
    if vtype == "ins":
        return INS_TARGETS
    if vtype == "invdel":
        return INVDEL_TARGETS
    return SNVINDEL_TARGETS


def ref_for(vtype: str) -> str:
    """Return the reference path (relative to repo root) for a variant type."""
    if vtype in ("sv", "ins", "invdel"):
        return "docs/benchmarking/sv/data/ref.fa"
    return "docs/benchmarking/snvindel/data/ref.fa"


# ─── Sample specification ─────────────────────────────────────────────────────

def make_sample(
    name: str,
    vtype: str,
    vaf: float,
    rep: str,
    coverage: float | None = None,
    family_size_mean: float | None = None,
    pcr_cycles: int | None = None,
    fragment_mean: float | None = None,
    seed: int | None = None,
) -> dict[str, Any]:
    """Build a sample specification dict."""
    tag = vaf_tag(vaf)
    cov    = coverage         if coverage         is not None else DEFAULTS["coverage"]
    fam    = family_size_mean if family_size_mean is not None else DEFAULTS["family_size_mean"]
    pcr    = pcr_cycles       if pcr_cycles       is not None else DEFAULTS["pcr_cycles"]
    fmean  = fragment_mean    if fragment_mean    is not None else DEFAULTS["fragment_mean"]
    s      = seed             if seed             is not None else seed_for(vtype, rep, vaf)

    if vtype in ("sv", "ins", "invdel"):
        truth_vcf = truth_vcf_sv(vtype, tag, rep)
    else:
        truth_vcf = truth_vcf_snvindel(vtype, tag, rep)

    return {
        "name":             name,
        "vtype":            vtype,
        "vaf":              vaf,
        "tag":              tag,
        "rep":              rep,
        "coverage":         cov,
        "family_size_mean": fam,
        "family_size_sd":   DEFAULTS["family_size_sd"],
        "pcr_cycles":       pcr,
        "fragment_mean":    fmean,
        "fragment_sd":      DEFAULTS["fragment_sd"],
        "seed":             s,
        "truth_vcf":        truth_vcf,
        "ref":              ref_for(vtype),
        "targets":          targets_for(vtype),
    }


# ─── Batch definitions ────────────────────────────────────────────────────────

def build_sample_list() -> list[dict[str, Any]]:
    """Return all samples to generate across all batches."""
    samples: list[dict[str, Any]] = []

    # Batch 1: new replicates c, d, e for SNV and indel
    for vtype in ("snv", "indel"):
        for rep in ("c", "d", "e"):
            for vaf in ALL_VAFS:
                tag = vaf_tag(vaf)
                name = f"{vtype}_vaf{tag}_{rep}"
                samples.append(make_sample(name, vtype, vaf, rep))

    # Batch 2: coverage sweep for SNV and indel
    for vtype in ("snv", "indel"):
        for cov in (1000, 2000, 8000, 12000):
            for vaf in SWEEP_VAFS:
                tag = vaf_tag(vaf)
                name = f"{vtype}_vaf{tag}_cov{int(cov)}_a"
                samples.append(make_sample(name, vtype, vaf, "a", coverage=float(cov)))

    # Batch 3: family size sweep for SNV and indel
    for vtype in ("snv", "indel"):
        for fam in (2.0, 6.0, 8.0):
            for vaf in SWEEP_VAFS:
                tag = vaf_tag(vaf)
                fam_str = f"fam{int(fam)}"
                name = f"{vtype}_vaf{tag}_{fam_str}_a"
                samples.append(make_sample(name, vtype, vaf, "a", family_size_mean=fam))

    # Batch 4: PCR cycles sweep for SNV and indel
    for vtype in ("snv", "indel"):
        for pcr in (5, 10, 12):
            for vaf in SWEEP_VAFS:
                tag = vaf_tag(vaf)
                name = f"{vtype}_vaf{tag}_pcr{pcr}_a"
                samples.append(make_sample(name, vtype, vaf, "a", pcr_cycles=pcr))

    # Batch 5: fragment size sweep for SNV and indel
    for vtype in ("snv", "indel"):
        for fmean in (140, 200, 250):
            for vaf in SWEEP_VAFS:
                tag = vaf_tag(vaf)
                name = f"{vtype}_vaf{tag}_frag{fmean}_a"
                samples.append(make_sample(name, vtype, vaf, "a", fragment_mean=float(fmean)))

    # Batch 6: SV new replicates c, d, e
    for rep in ("c", "d", "e"):
        for vaf in ALL_VAFS:
            tag = vaf_tag(vaf)
            name = f"sv_vaf{tag}_{rep}"
            samples.append(make_sample(name, "sv", vaf, rep))

    # Batch 7: INS new replicates c, d, e
    for rep in ("c", "d", "e"):
        for vaf in ALL_VAFS:
            tag = vaf_tag(vaf)
            name = f"ins_vaf{tag}_{rep}"
            samples.append(make_sample(name, "ins", vaf, rep))

    # Batch 8: INVDEL new replicates c, d, e
    for rep in ("c", "d", "e"):
        for vaf in ALL_VAFS:
            tag = vaf_tag(vaf)
            name = f"invdel_vaf{tag}_{rep}"
            samples.append(make_sample(name, "invdel", vaf, rep))

    return samples


# ─── Config generation ────────────────────────────────────────────────────────

def write_config(s: dict[str, Any]) -> Path:
    """Write a varforge YAML config for the sample. Returns the config path."""
    CONFIGS_DIR.mkdir(parents=True, exist_ok=True)
    config_path = CONFIGS_DIR / f"{s['name']}.yaml"

    sim_output_dir = f"docs/benchmarking/ml/results/sim_{s['name']}"
    sample_label   = s["name"].upper()
    vaf_pct        = s["vaf"] * 100.0

    content = (
        f"# varforge config: {sample_label} — VAF {vaf_pct:.3f}%\n"
        f"reference: {s['ref']}\n"
        f"\n"
        f"output:\n"
        f"  directory: {sim_output_dir}\n"
        f"  fastq: true\n"
        f"  bam: false\n"
        f"  truth_vcf: true\n"
        f"  manifest: true\n"
        f"\n"
        f"sample:\n"
        f"  name: {sample_label}\n"
        f"  read_length: 150\n"
        f"  coverage: {s['coverage']:.1f}\n"
        f"  platform: illumina\n"
        f"\n"
        f"fragment:\n"
        f"  model: normal\n"
        f"  mean: {s['fragment_mean']:.1f}\n"
        f"  sd: {s['fragment_sd']:.1f}\n"
        f"\n"
        f"quality:\n"
        f"  mean_quality: 37\n"
        f"  tail_decay: 0.002\n"
        f"\n"
        f"tumour:\n"
        f"  purity: {s['vaf'] * 2.0:.6f}\n"
        f"  ploidy: 2\n"
        f"\n"
        f"mutations:\n"
        f"  vcf: {s['truth_vcf'].relative_to(REPO)}\n"
        f"\n"
        f"umi:\n"
        f"  length: 5\n"
        f"  duplex: true\n"
        f"  pcr_cycles: {s['pcr_cycles']}\n"
        f"  family_size_mean: {s['family_size_mean']:.1f}\n"
        f"  family_size_sd: {s['family_size_sd']:.1f}\n"
        f"\n"
        f"chromosomes:\n"
        f"  - chr1\n"
        f"\n"
        f"seed: {s['seed']}\n"
    )
    config_path.write_text(content)
    return config_path


def write_varforge_cmd(s: dict[str, Any], config_path: Path) -> Path:
    """Write the varforge command file. Returns the path."""
    cmd_path = CONFIGS_DIR / f"{s['name']}_cmd.txt"
    config_rel = config_path.relative_to(REPO)
    cmd_path.write_text(
        f"# Run from the repository root.\n"
        f"varforge simulate --config {config_rel}\n"
    )
    return cmd_path


# ─── Simulation ───────────────────────────────────────────────────────────────

def run_simulation(s: dict[str, Any], config_path: Path) -> bool:
    """Run varforge simulate. Returns True on success."""
    sim_dir = RESULTS_DIR / f"sim_{s['name']}"
    # Check if already done — varforge writes R1/R2 FASTQs.
    if sim_dir.exists() and any(sim_dir.glob("*_R1.fastq.gz")):
        return True

    sim_dir.mkdir(parents=True, exist_ok=True)
    config_rel = str(config_path.relative_to(REPO))

    print(f"[SIM]  simulating {s['name']}", flush=True)
    result = subprocess.run(
        ["varforge", "simulate", "--config", config_rel],
        cwd=str(REPO),
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        print(f"[WARN] varforge failed for {s['name']}: {result.stderr[-500:]}", flush=True)
        return False

    # Decompress truth VCF — kam cannot read .vcf.gz directly.
    for gz_vcf in sim_dir.glob("*.truth.vcf.gz"):
        plain_vcf = gz_vcf.with_suffix("")  # strips .gz -> .truth.vcf
        if not plain_vcf.exists():
            subprocess.run(
                ["gunzip", "-k", str(gz_vcf)],
                cwd=str(REPO),
                check=False,
            )

    return True


# ─── kam pipeline ─────────────────────────────────────────────────────────────

def run_kam(s: dict[str, Any]) -> bool:
    """Run kam discovery and tumour-informed on a simulated sample.

    Returns True if both modes completed.
    """
    sim_dir = RESULTS_DIR / f"sim_{s['name']}"
    out_dir = RESULTS_DIR / f"kam_{s['name']}"

    disc_vcf = out_dir / "calls_discovery.vcf"
    ti_vcf   = out_dir / "calls_tumour_informed.vcf"

    if disc_vcf.exists() and ti_vcf.exists():
        return True

    # Locate FASTQs.
    r1_files = sorted(sim_dir.glob("*_R1.fastq.gz"))
    r2_files = sorted(sim_dir.glob("*_R2.fastq.gz"))
    if not r1_files or not r2_files:
        print(f"[WARN] No FASTQs in {sim_dir}", flush=True)
        return False

    r1 = str(r1_files[0])
    r2 = str(r2_files[0])

    # Locate varforge truth VCF (uncompressed — .vcf.gz was decompressed after sim).
    truth_files = sorted(sim_dir.glob("*.truth.vcf"))
    if not truth_files:
        # Fall back to gzipped if decompression hasn't happened yet.
        truth_files_gz = sorted(sim_dir.glob("*.truth.vcf.gz"))
        if truth_files_gz:
            gz = truth_files_gz[0]
            plain = gz.with_suffix("")
            subprocess.run(["gunzip", "-k", str(gz)], cwd=str(REPO), check=False)
            truth_files = [plain] if plain.exists() else []
    if not truth_files:
        print(f"[WARN] No truth VCF in {sim_dir}", flush=True)
        return False
    truth_vcf = str(truth_files[0])

    targets = str(s["targets"])
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"[KAM]  discovery  {s['name']}", flush=True)

    # SV types use position-tolerance for tumour-informed matching.
    is_sv = s["vtype"] in ("sv", "ins", "invdel")

    tmp_disc = out_dir / "tmp_disc"
    disc_cmd = [
        str(KAM), "run",
        "--r1", r1, "--r2", r2,
        "--targets", targets,
        "--output-dir", str(tmp_disc),
        "--output-format-override", "vcf,tsv",
    ]
    result = subprocess.run(disc_cmd, cwd=str(REPO), capture_output=True, text=True)
    if result.returncode != 0:
        print(f"[WARN] kam discovery failed for {s['name']}: {result.stderr[-300:]}", flush=True)
        return False
    shutil.copy2(tmp_disc / "variants.vcf", disc_vcf)
    shutil.copy2(tmp_disc / "variants.tsv", out_dir / "calls_discovery.tsv")
    shutil.rmtree(tmp_disc, ignore_errors=True)

    print(f"[KAM]  ti         {s['name']}", flush=True)
    tmp_ti = out_dir / "tmp_ti"
    ti_cmd = [
        str(KAM), "run",
        "--r1", r1, "--r2", r2,
        "--targets", targets,
        "--target-variants", truth_vcf,
        "--output-dir", str(tmp_ti),
        "--output-format-override", "vcf,tsv",
    ]
    if is_sv:
        ti_cmd += ["--ti-position-tolerance", "10"]
    result = subprocess.run(ti_cmd, cwd=str(REPO), capture_output=True, text=True)
    if result.returncode != 0:
        print(f"[WARN] kam tumour-informed failed for {s['name']}: {result.stderr[-300:]}", flush=True)
        return False
    shutil.copy2(tmp_ti / "variants.vcf", ti_vcf)
    shutil.copy2(tmp_ti / "variants.tsv", out_dir / "calls_tumour_informed.tsv")
    shutil.rmtree(tmp_ti, ignore_errors=True)

    return True


# ─── VCF parsing ──────────────────────────────────────────────────────────────

def parse_info(info_str: str) -> dict[str, str]:
    """Parse a VCF INFO field into a dict."""
    result: dict[str, str] = {}
    for field in info_str.split(";"):
        if "=" in field:
            k, v = field.split("=", 1)
            result[k] = v
        else:
            result[field] = ""
    return result


def open_vcf(vcf_path: Path):
    """Open a VCF file, handling both plain text and gzip."""
    if vcf_path.suffix == ".gz":
        return gzip.open(vcf_path, "rt")
    return vcf_path.open()


def vcf_to_truth_tsv(vcf_path: Path, out_path: Path) -> None:
    """Convert a truth VCF to a simple TSV.

    Columns: chrom, pos, ref, alt, vaf, type
    """
    rows = ["chrom\tpos\tref\talt\tvaf\ttype"]
    with open_vcf(vcf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            chrom, pos, _, ref, alt = parts[0], parts[1], parts[2], parts[3], parts[4]
            info = parse_info(parts[7]) if len(parts) > 7 else {}
            vaf          = info.get("VAF", "")
            variant_type = info.get("SVTYPE", info.get("TYPE", ""))
            rows.append(f"{chrom}\t{pos}\t{ref}\t{alt}\t{vaf}\t{variant_type}")
    out_path.write_text("\n".join(rows) + "\n")


def vcf_to_calls_tsv(vcf_path: Path, out_path: Path) -> None:
    """Convert a kam calls VCF to a flat TSV.

    Columns: chrom, pos, ref, alt, filter,
             vaf, vaf_lo, vaf_hi, nref, nalt, ndupalt, nsimalt, sbp, conf
    """
    header = (
        "chrom\tpos\tref\talt\tfilter\t"
        "vaf\tvaf_lo\tvaf_hi\tnref\tnalt\tndupalt\tnsimalt\tsbp\tconf"
    )
    rows = [header]
    with vcf_path.open() as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue
            chrom, pos = parts[0], parts[1]
            ref, alt   = parts[3], parts[4]
            filt       = parts[6]
            info = parse_info(parts[7])
            vaf     = info.get("VAF", "")
            vaf_lo  = info.get("VAF_LO", "")
            vaf_hi  = info.get("VAF_HI", "")
            nref    = info.get("NREF", "")
            nalt    = info.get("NALT", "")
            ndupalt = info.get("NDUPALT", "")
            nsimalt = info.get("NSIMALT", "")
            sbp     = info.get("SBP", "")
            conf    = info.get("CONF", "")
            rows.append(
                f"{chrom}\t{pos}\t{ref}\t{alt}\t{filt}\t"
                f"{vaf}\t{vaf_lo}\t{vaf_hi}\t{nref}\t{nalt}\t"
                f"{ndupalt}\t{nsimalt}\t{sbp}\t{conf}"
            )
    out_path.write_text("\n".join(rows) + "\n")


# ─── Per-sample directory ─────────────────────────────────────────────────────

def build_sample_dir(s: dict[str, Any], config_path: Path, cmd_path: Path) -> bool:
    """Build the per-sample directory. Returns True on success."""
    sample_dir = SAMPLES_DIR / s["name"]

    # Skip if already fully built.
    if (sample_dir / "discovery.tsv").exists():
        return True

    out_dir = RESULTS_DIR / f"kam_{s['name']}"
    disc_vcf = out_dir / "calls_discovery.vcf"
    ti_vcf   = out_dir / "calls_tumour_informed.vcf"
    sim_dir  = RESULTS_DIR / f"sim_{s['name']}"

    # Locate varforge-written truth VCF (uncompressed).
    truth_files = sorted(sim_dir.glob("*.truth.vcf"))
    if not truth_files:
        # Try gzipped as fallback and decompress.
        for gz in sorted(sim_dir.glob("*.truth.vcf.gz")):
            plain = gz.with_suffix("")
            if not plain.exists():
                plain.write_bytes(gzip.open(gz, "rb").read())
            truth_files.append(plain)

    if not truth_files or not disc_vcf.exists() or not ti_vcf.exists():
        return False

    varforge_truth_vcf = truth_files[0]

    sample_dir.mkdir(parents=True, exist_ok=True)

    # config.yaml
    shutil.copy2(config_path, sample_dir / "config.yaml")

    # varforge_cmd.txt
    config_rel = config_path.relative_to(REPO)
    (sample_dir / "varforge_cmd.txt").write_text(
        f"# Run from the repository root.\n"
        f"varforge simulate --config {config_rel}\n"
    )

    # truth.tsv from varforge output truth VCF
    vcf_to_truth_tsv(varforge_truth_vcf, sample_dir / "truth.tsv")

    # discovery.tsv and tumour_informed.tsv
    vcf_to_calls_tsv(disc_vcf, sample_dir / "discovery.tsv")
    vcf_to_calls_tsv(ti_vcf,   sample_dir / "tumour_informed.tsv")

    # params.json
    params = {
        "coverage":         s["coverage"],
        "family_size_mean": s["family_size_mean"],
        "pcr_cycles":       s["pcr_cycles"],
        "fragment_mean":    s["fragment_mean"],
        "vaf":              s["vaf"],
        "replicate":        s["rep"],
        "variant_type":     s["vtype"],
    }
    (sample_dir / "params.json").write_text(json.dumps(params, indent=2) + "\n")

    return True


# ─── Process one sample ───────────────────────────────────────────────────────

def process_sample(s: dict[str, Any]) -> tuple[str, str]:
    """Run the full pipeline for one sample. Returns (name, status)."""
    name = s["name"]

    # Skip if already fully complete.
    if (SAMPLES_DIR / name / "discovery.tsv").exists():
        return name, "skipped"

    # Verify truth VCF exists before doing anything.
    if not s["truth_vcf"].exists():
        return name, f"missing_truth:{s['truth_vcf']}"

    try:
        config_path = write_config(s)
        cmd_path    = write_varforge_cmd(s, config_path)

        if not run_simulation(s, config_path):
            return name, "sim_failed"

        if not run_kam(s):
            return name, "kam_failed"

        if not build_sample_dir(s, config_path, cmd_path):
            return name, "build_failed"

        print(f"[DONE] {name}", flush=True)
        return name, "done"

    except Exception as exc:  # noqa: BLE001
        return name, f"error:{exc}"


# ─── Main ─────────────────────────────────────────────────────────────────────

def main() -> None:
    if not KAM.exists():
        print(f"[ERROR] kam binary not found: {KAM}", file=sys.stderr)
        sys.exit(1)

    # Create output directories.
    for d in (CONFIGS_DIR, RESULTS_DIR, SAMPLES_DIR):
        d.mkdir(parents=True, exist_ok=True)

    samples = build_sample_list()
    print(f"Total samples to process: {len(samples)}", flush=True)

    counts: dict[str, int] = {}

    with ThreadPoolExecutor(max_workers=4) as pool:
        futures = {pool.submit(process_sample, s): s["name"] for s in samples}
        for fut in as_completed(futures):
            name, status = fut.result()
            key = status if status in ("done", "skipped") else "failed"
            counts[key] = counts.get(key, 0) + 1
            if key == "failed":
                print(f"[FAIL] {name}: {status}", flush=True)

    print()
    print("=== SUMMARY ===")
    print(f"  Done:    {counts.get('done', 0)}")
    print(f"  Skipped: {counts.get('skipped', 0)}")
    print(f"  Failed:  {counts.get('failed', 0)}")
    total_complete = counts.get("done", 0) + counts.get("skipped", 0)
    print(f"  Total complete: {total_complete} / {len(samples)}")


if __name__ == "__main__":
    main()
