#!/usr/bin/env python3
"""Pipeline runner for the Twist duplex ML dataset.

For each of 11,000 samples (10k train + 1k test):
  1. Run varforge simulate → generates FASTQs and a truth VCF.
  2. Run kam discovery mode → calls_discovery.tsv / .vcf.
  3. Run kam tumour-informed mode → calls_tumour_informed.tsv / .vcf.
  4. Write params.json to the sample directory.

Per-sample output:
  docs/benchmarking/ml-twist-duplex/simulations/{split}/{sample_name}/
    {sample_name}_R1.fastq.gz, {sample_name}_R2.fastq.gz
    {sample_name}.truth.vcf.gz      (from varforge)
    calls_discovery.tsv / .vcf      (from kam)
    calls_tumour_informed.tsv / .vcf (from kam)
    params.json

Processes in batches of 500. After each batch, tars all sample directories
and uploads to Nextcloud, then continues. Checkpointing allows resuming.

Usage:
    python3 scripts/ml/run_twist_duplex_pipeline.py [--split train|test|all]
                                                     [--workers N]
                                                     [--batch-size N]
                                                     [--dry-run]
                                                     [--skip-upload]

Run from the repository root.
"""

from __future__ import annotations

import argparse
import json
import logging
import subprocess
import sys
import tarfile
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

# ─── Paths ────────────────────────────────────────────────────────────────────

REPO = Path(__file__).resolve().parent.parent.parent
ML_DIR = REPO / "docs" / "benchmarking" / "ml-twist-duplex"
MANIFEST = ML_DIR / "manifest.json"
CHECKPOINT = ML_DIR / "checkpoint.json"
SIM_DIR = ML_DIR / "simulations"
KAM = REPO / "target" / "release" / "kam"

# ─── Targets FAs ─────────────────────────────────────────────────────────────

SNVINDEL_TARGETS = REPO / "docs/benchmarking/snvindel/data/snvindel_targets.fa"
SV_TARGETS = REPO / "docs/benchmarking/sv/data/sv_suite_targets.fa"
INS_TARGETS = REPO / "docs/benchmarking/sv/data/ins_targets.fa"
INVDEL_TARGETS = REPO / "docs/benchmarking/sv/data/invdel_targets.fa"

# SV types needing position tolerance for tumour-informed matching.
SV_TYPES = {"sv_dupinv", "ins", "invdel"}

# ─── Nextcloud ────────────────────────────────────────────────────────────────

NEXTCLOUD_URL = "https://nextcloudlocal.trentz.me/public.php/dav/files/pTizAiSAJQsPcDo"
NEXTCLOUD_TOKEN = "pTizAiSAJQsPcDo"

# ─── Parallelism ─────────────────────────────────────────────────────────────

# 8 workers × 2 threads each = 16 threads total, matching a 16-thread machine.
# More workers (less threads per process) is better here: the reference is small
# (2 kbp) so per-process parallelism gives diminishing returns, while more
# concurrent samples saturates the CPU more evenly.
DEFAULT_WORKERS = 8
PROCESS_THREADS = 2
DEFAULT_BATCH_SIZE = 500

# ─── Logging ─────────────────────────────────────────────────────────────────

logging.basicConfig(
    level=logging.INFO,
    format="[%(levelname)s] %(message)s",
    stream=sys.stderr,
)
log = logging.getLogger(__name__)


# ─── Helpers ─────────────────────────────────────────────────────────────────

def targets_for(vtype: str) -> Path:
    """Return the targets FASTA for a variant type.

    Args:
        vtype: Variant type key from the manifest.

    Returns:
        Path to the targets FASTA file.
    """
    if vtype in ("snv", "indel"):
        return SNVINDEL_TARGETS
    if vtype == "sv_dupinv":
        return SV_TARGETS
    if vtype == "ins":
        return INS_TARGETS
    if vtype == "invdel":
        return INVDEL_TARGETS
    return SNVINDEL_TARGETS


def find_fastqs(sim_dir: Path, name: str) -> tuple[Path | None, Path | None]:
    """Locate R1 and R2 FASTQ files in a varforge output directory.

    Args:
        sim_dir: Varforge output directory for this sample.
        name: Sample name.

    Returns:
        (r1, r2) paths, or (None, None) if not found.
    """
    r1 = sim_dir / f"{name}_R1.fastq.gz"
    r2 = sim_dir / f"{name}_R2.fastq.gz"
    if r1.exists() and r2.exists():
        return r1, r2
    # Fallback: glob search.
    r1_matches = sorted(sim_dir.glob("*_R1*.fastq.gz"))
    r2_matches = sorted(sim_dir.glob("*_R2*.fastq.gz"))
    if r1_matches and r2_matches:
        return r1_matches[0], r2_matches[0]
    return None, None


# ─── Pipeline steps ───────────────────────────────────────────────────────────

def run_varforge(config: str, name: str) -> bool:
    """Run varforge simulate for one sample.

    Args:
        config: Relative path to the varforge YAML config.
        name: Sample name (for logging).

    Returns:
        True on success, False on failure.
    """
    result = subprocess.run(
        ["varforge", "simulate", "--config", config, "-t", str(PROCESS_THREADS)],
        cwd=str(REPO),
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        log.error("varforge failed for %s (rc=%d): %s",
                  name, result.returncode, result.stderr[-300:])
        return False
    return True


def run_kam_discovery(
    name: str,
    r1: Path,
    r2: Path,
    targets: Path,
    out_dir: Path,
) -> bool:
    """Run kam in discovery mode.

    Args:
        name: Sample name (for logging).
        r1: R1 FASTQ path.
        r2: R2 FASTQ path.
        targets: Targets FASTA path.
        out_dir: Sample output directory.

    Returns:
        True on success, False on failure.
    """
    tmp = out_dir / "tmp_disc"
    tmp.mkdir(parents=True, exist_ok=True)

    result = subprocess.run(
        [
            str(KAM), "run",
            "--r1", str(r1),
            "--r2", str(r2),
            "--targets", str(targets),
            "--output-dir", str(tmp),
            "--output-format-override", "vcf,tsv",
            "--threads", str(PROCESS_THREADS),
        ],
        cwd=str(REPO),
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        log.error("kam discovery failed for %s (rc=%d): %s",
                  name, result.returncode, result.stderr[-300:])
        return False

    tsv = tmp / "variants.tsv"
    vcf = tmp / "variants.vcf"
    if not tsv.exists():
        log.error("kam discovery produced no TSV for %s", name)
        return False

    tsv.rename(out_dir / "calls_discovery.tsv")
    if vcf.exists():
        vcf.rename(out_dir / "calls_discovery.vcf")
    # Remove tmp dir (may still contain stray files; ignore errors).
    try:
        for f in tmp.iterdir():
            f.unlink(missing_ok=True)
        tmp.rmdir()
    except OSError:
        pass
    return True


def run_kam_tumour_informed(
    name: str,
    r1: Path,
    r2: Path,
    targets: Path,
    truth_vcf: Path,
    vtype: str,
    out_dir: Path,
) -> bool:
    """Run kam in tumour-informed mode.

    Args:
        name: Sample name (for logging).
        r1: R1 FASTQ path.
        r2: R2 FASTQ path.
        targets: Targets FASTA path.
        truth_vcf: Truth VCF for target-variants filter.
        vtype: Variant type key (SV types get position tolerance).
        out_dir: Sample output directory.

    Returns:
        True on success, False on failure.
    """
    tmp = out_dir / "tmp_ti"
    tmp.mkdir(parents=True, exist_ok=True)

    cmd = [
        str(KAM), "run",
        "--r1", str(r1),
        "--r2", str(r2),
        "--targets", str(targets),
        "--output-dir", str(tmp),
        "--output-format-override", "vcf,tsv",
        "--target-variants", str(truth_vcf),
        "--threads", str(PROCESS_THREADS),
    ]
    if vtype in SV_TYPES:
        cmd += ["--ti-position-tolerance-override", "10"]

    result = subprocess.run(
        cmd,
        cwd=str(REPO),
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        log.error("kam tumour-informed failed for %s (rc=%d): %s",
                  name, result.returncode, result.stderr[-300:])
        return False

    tsv = tmp / "variants.tsv"
    vcf = tmp / "variants.vcf"
    if not tsv.exists():
        log.error("kam tumour-informed produced no TSV for %s", name)
        return False

    tsv.rename(out_dir / "calls_tumour_informed.tsv")
    if vcf.exists():
        vcf.rename(out_dir / "calls_tumour_informed.vcf")
    try:
        for f in tmp.iterdir():
            f.unlink(missing_ok=True)
        tmp.rmdir()
    except OSError:
        pass
    return True


def write_params(sample: dict, out_dir: Path) -> None:
    """Write params.json for feature extraction.

    Args:
        sample: Manifest entry for this sample.
        out_dir: Sample output directory.
    """
    params = {
        "name": sample["name"],
        "split": sample["split"],
        "vtype": sample["vtype"],
        "coverage": sample["coverage"],
        "family_size_mean": sample["family_size_mean"],
        "family_size_sd": sample["family_size_sd"],
        "pcr_cycles": sample["pcr_cycles"],
        "fragment_mean": sample["fragment_mean"],
        "fragment_sd": sample["fragment_sd"],
        "mean_quality": sample["mean_quality"],
        "vaf": sample["vaf"],
        "purity": sample["purity"],
        "seed": sample["seed"],
        "truth_vcf": sample["truth_vcf"],
    }
    (out_dir / "params.json").write_text(json.dumps(params, indent=2))


# ─── Per-sample pipeline ──────────────────────────────────────────────────────

def process_sample(sample: dict) -> tuple[str, str]:
    """Run the full pipeline for one sample.

    Args:
        sample: Manifest entry for this sample.

    Returns:
        (sample_name, status) where status is 'ok' or an error description.
    """
    name = sample["name"]
    vtype = sample["vtype"]
    config = sample["config"]
    truth_vcf = REPO / sample["truth_vcf"]
    out_dir = REPO / sample["sim_dir"]
    out_dir.mkdir(parents=True, exist_ok=True)

    targets = targets_for(vtype)

    # Step 1: varforge simulation.
    if not run_varforge(config, name):
        return name, "varforge_failed"

    # Step 2: locate FASTQs.
    r1, r2 = find_fastqs(out_dir, name)
    if r1 is None:
        log.error("No FASTQs found for %s in %s", name, out_dir)
        return name, "no_fastqs"

    # Step 3: kam discovery.
    if not run_kam_discovery(name, r1, r2, targets, out_dir):
        return name, "discovery_failed"

    # Step 4: kam tumour-informed.
    if not run_kam_tumour_informed(name, r1, r2, targets, truth_vcf, vtype, out_dir):
        return name, "tumour_informed_failed"

    # Step 5: write params.json.
    write_params(sample, out_dir)

    return name, "ok"


# ─── Checkpoint ───────────────────────────────────────────────────────────────

def load_checkpoint() -> set[str]:
    """Return set of already-completed sample names."""
    if CHECKPOINT.exists():
        data = json.loads(CHECKPOINT.read_text())
        return set(data.get("completed", []))
    return set()


def save_checkpoint(completed: set[str]) -> None:
    """Persist completed sample names to disk."""
    CHECKPOINT.write_text(json.dumps({"completed": sorted(completed)}, indent=2))


# ─── Nextcloud upload ─────────────────────────────────────────────────────────

def upload_batch(
    split: str,
    batch_idx: int,
    sample_names: list[str],
    dry_run: bool = False,
) -> None:
    """Tar and upload a batch of sample directories to Nextcloud.

    Args:
        split: 'train' or 'test'.
        batch_idx: Batch number (0-indexed, global).
        sample_names: Names of completed samples in this batch.
        dry_run: If True, skip actual upload.
    """
    split_dir = SIM_DIR / split
    tarball = Path(f"/tmp/ml_twist_{split}_batch_{batch_idx:03d}.tar.gz")
    nc_path = f"benchmarking/ml-twist-duplex/{split}/batch_{batch_idx:03d}.tar.gz"

    log.info("Tarballing batch %d (%d samples)...", batch_idx, len(sample_names))
    with tarfile.open(tarball, "w:gz") as tf:
        for sname in sample_names:
            d = split_dir / sname
            if d.exists():
                tf.add(d, arcname=sname)

    size_mb = tarball.stat().st_size / 1_048_576
    log.info("Tarball %.1f MB — uploading...", size_mb)

    if dry_run:
        log.info("[DRY RUN] Would upload to %s/%s", NEXTCLOUD_URL, nc_path)
        tarball.unlink(missing_ok=True)
        return

    result = subprocess.run(
        [
            "curl", "-T", str(tarball),
            f"{NEXTCLOUD_URL}/{nc_path}",
            "-u", f"{NEXTCLOUD_TOKEN}:",
            "--silent", "--show-error",
        ],
        capture_output=True,
        text=True,
    )
    tarball.unlink(missing_ok=True)

    if result.returncode != 0:
        log.error("Upload failed for batch %d: %s", batch_idx, result.stderr)
    else:
        log.info("Batch %d uploaded.", batch_idx)


# ─── Split processing ─────────────────────────────────────────────────────────

def process_split(
    samples: list[dict],
    split: str,
    workers: int,
    batch_size: int,
    dry_run: bool,
    skip_upload: bool,
) -> None:
    """Process all samples for one split.

    Args:
        samples: Manifest entries for this split.
        split: 'train' or 'test'.
        workers: Concurrent processes.
        batch_size: Samples per upload batch.
        dry_run: Log only, no varforge/kam/upload.
        skip_upload: Run varforge/kam but skip upload.
    """
    completed = load_checkpoint()
    remaining = [s for s in samples if s["name"] not in completed]
    n_already = len(samples) - len(remaining)

    log.info("Total samples: %d | Already done: %d | To run: %d",
             len(samples), n_already, len(remaining))

    batches = [remaining[i: i + batch_size] for i in range(0, len(remaining), batch_size)]
    log.info("Batches: %d × ~%d samples", len(batches), batch_size)

    for b_local, batch in enumerate(batches):
        b_global = n_already // batch_size + b_local
        log.info("[BATCH %d] Processing %d samples...", b_global, len(batch))
        t0 = time.monotonic()
        ok: list[str] = []
        failed: list[str] = []

        if dry_run:
            for s in batch:
                log.info("  [DRY RUN] %s", s["name"])
            ok = [s["name"] for s in batch]
        else:
            with ThreadPoolExecutor(max_workers=workers) as pool:
                futures = {pool.submit(process_sample, s): s["name"] for s in batch}
                for fut in as_completed(futures):
                    name, status = fut.result()
                    if status == "ok":
                        ok.append(name)
                        completed.add(name)
                        save_checkpoint(completed)
                    else:
                        failed.append(name)
                        log.warning("FAILED %s: %s", name, status)

        elapsed = time.monotonic() - t0
        log.info("[BATCH %d] %.0fs — %d ok, %d failed",
                 b_global, elapsed, len(ok), len(failed))
        if failed:
            log.warning("Failed: %s", failed)

        if ok and not skip_upload:
            upload_batch(split, b_global, ok, dry_run=dry_run)

    log.info("Split '%s' complete.", split)


# ─── Entry point ─────────────────────────────────────────────────────────────

def main() -> None:
    """Parse arguments and run the pipeline."""
    parser = argparse.ArgumentParser(
        description="Run varforge + kam for all Twist duplex ML dataset samples.",
    )
    parser.add_argument(
        "--split", choices=["train", "test", "all"], default="all",
        help="Which split to process (default: all).",
    )
    parser.add_argument(
        "--workers", type=int, default=DEFAULT_WORKERS,
        help=f"Concurrent sample processes (default: {DEFAULT_WORKERS}).",
    )
    parser.add_argument(
        "--batch-size", type=int, default=DEFAULT_BATCH_SIZE,
        help=f"Samples per upload batch (default: {DEFAULT_BATCH_SIZE}).",
    )
    parser.add_argument(
        "--dry-run", action="store_true",
        help="Log without running anything.",
    )
    parser.add_argument(
        "--skip-upload", action="store_true",
        help="Run pipeline but skip Nextcloud upload.",
    )
    args = parser.parse_args()

    if not MANIFEST.exists():
        log.error("Manifest not found: %s — run generate_twist_duplex_configs.py first.", MANIFEST)
        sys.exit(1)

    if not KAM.exists():
        log.error("kam binary not found: %s — run `cargo build --release` first.", KAM)
        sys.exit(1)

    manifest = json.loads(MANIFEST.read_text())
    all_samples = manifest["samples"]

    SIM_DIR.mkdir(parents=True, exist_ok=True)
    splits: list[tuple[str, list[dict]]] = []
    if args.split in ("train", "all"):
        splits.append(("train", [s for s in all_samples if s["split"] == "train"]))
    if args.split in ("test", "all"):
        splits.append(("test", [s for s in all_samples if s["split"] == "test"]))

    for split, samples in splits:
        log.info("=== Split: %s ===", split)
        (SIM_DIR / split).mkdir(parents=True, exist_ok=True)
        process_split(
            samples=samples,
            split=split,
            workers=args.workers,
            batch_size=args.batch_size,
            dry_run=args.dry_run,
            skip_upload=args.skip_upload,
        )

    log.info("Done.")


if __name__ == "__main__":
    main()
