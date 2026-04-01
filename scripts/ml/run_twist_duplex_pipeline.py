#!/usr/bin/env python3
"""Full pipeline runner for the Twist duplex ML dataset.

For each of 11,000 samples (10k train + 1k test):
  1. Run varforge simulate to generate FASTQs and truth VCF.
  2. Run kam (discovery mode) on the FASTQs.
  3. Run kam (tumour-informed mode) on the FASTQs.
  4. Write params.json to the sample directory.
  5. Copy the truth VCF to the sample directory.

Processes in batches of 500. After each batch, uploads to Nextcloud and
logs progress. Checkpointing allows resuming interrupted runs.

Usage:
    python3 scripts/ml/run_twist_duplex_pipeline.py [--split train|test|all]
                                                      [--workers N]
                                                      [--batch-size N]
                                                      [--dry-run]

Run from the repository root.
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any

REPO = Path(__file__).resolve().parent.parent.parent

# ─── Paths ────────────────────────────────────────────────────────────────────

ML_DIR = REPO / "docs" / "benchmarking" / "ml-twist-duplex"
CONFIGS_TRAIN = ML_DIR / "configs" / "train"
CONFIGS_TEST = ML_DIR / "configs" / "test"
SIM_DIR = ML_DIR / "simulations"
CHECKPOINT_PATH = ML_DIR / "checkpoint.json"
MANIFEST_PATH = ML_DIR / "manifest.json"

KAM = REPO / "target" / "release" / "kam"

# Targets FASTA: used for all variant types (the ML dataset uses a single
# synthetic chr1 that is short enough for all variant types).
# SNV/indel samples use the snvindel targets; SV samples use the SV targets.
SNVINDEL_TARGETS = REPO / "docs" / "benchmarking" / "snvindel" / "data" / "snvindel_targets.fa"
SV_TARGETS = REPO / "docs" / "benchmarking" / "sv" / "data" / "sv_suite_targets.fa"
INS_TARGETS = REPO / "docs" / "benchmarking" / "sv" / "data" / "ins_targets.fa"
INVDEL_TARGETS = REPO / "docs" / "benchmarking" / "sv" / "data" / "invdel_targets.fa"

# Nextcloud upload config
NEXTCLOUD_URL = "https://nextcloudlocal.trentz.me/public.php/dav/files/pTizAiSAJQsPcDo"
NEXTCLOUD_TOKEN = "pTizAiSAJQsPcDo"

DEFAULT_WORKERS = 12
DEFAULT_BATCH_SIZE = 500

SV_TYPES = {"sv_del", "sv_dup", "sv_inv", "sv_large_del", "ins", "invdel"}


# ─── Checkpoint ───────────────────────────────────────────────────────────────

def load_checkpoint() -> dict[str, str]:
    """Load the checkpoint file, returning a dict of sample_name -> status.

    Returns:
        Dict mapping sample name to 'done' or 'failed'.
    """
    if CHECKPOINT_PATH.exists():
        with open(CHECKPOINT_PATH) as fh:
            return json.load(fh)
    return {}


def save_checkpoint(checkpoint: dict[str, str]) -> None:
    """Write the checkpoint file atomically.

    Args:
        checkpoint: Dict of sample_name -> status.
    """
    tmp = CHECKPOINT_PATH.with_suffix(".tmp")
    with open(tmp, "w") as fh:
        json.dump(checkpoint, fh, indent=2)
    tmp.replace(CHECKPOINT_PATH)


# ─── Manifest loading ─────────────────────────────────────────────────────────

def load_manifest() -> list[dict[str, Any]]:
    """Load the dataset manifest, raising if not found.

    Returns:
        List of sample dicts from manifest.json.

    Raises:
        SystemExit: If manifest.json does not exist.
    """
    if not MANIFEST_PATH.exists():
        print(
            f"[ERROR] Manifest not found: {MANIFEST_PATH}\n"
            "Run scripts/ml/generate_twist_duplex_configs.py first.",
            file=sys.stderr,
        )
        sys.exit(1)
    with open(MANIFEST_PATH) as fh:
        data = json.load(fh)
    return data["samples"]


# ─── Target selection ─────────────────────────────────────────────────────────

def targets_for(vtype: str) -> Path:
    """Return the targets FASTA for a given variant type.

    Args:
        vtype: Variant type string, e.g. 'snv', 'ins'.

    Returns:
        Path to the targets FASTA file.
    """
    if vtype == "ins":
        return INS_TARGETS
    if vtype == "invdel":
        return INVDEL_TARGETS
    if vtype in SV_TYPES:
        return SV_TARGETS
    return SNVINDEL_TARGETS


# ─── FASTQ discovery ──────────────────────────────────────────────────────────

def find_fastqs(sample_dir: Path) -> tuple[Path | None, Path | None]:
    """Find R1 and R2 FASTQ files in a varforge output directory.

    Args:
        sample_dir: Directory to search.

    Returns:
        Tuple (r1_path, r2_path) or (None, None) if not found.
    """
    r1s = sorted(sample_dir.glob("*_R1.fastq.gz"))
    r2s = sorted(sample_dir.glob("*_R2.fastq.gz"))
    if r1s and r2s:
        return r1s[0], r2s[0]
    return None, None


def find_truth_vcf(sample_dir: Path) -> Path | None:
    """Find the truth VCF produced by varforge in a sample directory.

    Decompresses a .gz VCF if only the compressed form exists, as kam
    requires an uncompressed VCF for --target-variants.

    Args:
        sample_dir: varforge output directory.

    Returns:
        Path to the truth VCF, or None if not found.
    """
    import gzip as _gzip

    plain = sorted(sample_dir.glob("*.truth.vcf"))
    if plain:
        return plain[0]

    gz = sorted(sample_dir.glob("*.truth.vcf.gz"))
    if gz:
        out = gz[0].with_suffix("")  # strips .gz
        if not out.exists():
            with _gzip.open(gz[0], "rt") as fin, open(out, "w") as fout:
                fout.write(fin.read())
        return out

    return None


# ─── Pipeline steps ───────────────────────────────────────────────────────────

def run_varforge(config_path: Path, sample_name: str) -> bool:
    """Run varforge simulate for a single sample.

    Args:
        config_path: Path to the varforge YAML config.
        sample_name: Sample identifier (for logging).

    Returns:
        True on success, False on failure.
    """
    result = subprocess.run(
        ["varforge", "simulate", "--config", str(config_path.relative_to(REPO))],
        cwd=str(REPO),
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        print(
            f"[VARFORGE FAIL] {sample_name}: {result.stderr[-600:]}",
            file=sys.stderr,
            flush=True,
        )
        return False
    return True


def run_kam_discovery(
    sample_name: str,
    r1: Path,
    r2: Path,
    targets: Path,
    out_dir: Path,
) -> bool:
    """Run kam in discovery mode.

    Args:
        sample_name: For logging.
        r1: R1 FASTQ path.
        r2: R2 FASTQ path.
        targets: Targets FASTA path.
        out_dir: Directory where outputs will be written.

    Returns:
        True on success.
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
            "--output-format", "vcf,tsv",
        ],
        cwd=str(REPO),
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        print(
            f"[KAM DISC FAIL] {sample_name}: {result.stderr[-600:]}",
            file=sys.stderr,
            flush=True,
        )
        shutil.rmtree(tmp, ignore_errors=True)
        return False

    for src, dst_name in [
        (tmp / "variants.tsv", "discovery.tsv"),
        (tmp / "variants.vcf", "discovery.vcf"),
    ]:
        if src.exists():
            try:
                shutil.copy2(src, out_dir / dst_name)
            except OSError as exc:
                print(f"[WARN] copy failed {src}: {exc}", file=sys.stderr)

    shutil.rmtree(tmp, ignore_errors=True)
    return True


def run_kam_tumour_informed(
    sample_name: str,
    r1: Path,
    r2: Path,
    targets: Path,
    truth_vcf: Path,
    out_dir: Path,
    vtype: str,
) -> bool:
    """Run kam in tumour-informed mode.

    For SV types, passes --ti-position-tolerance 10 because the called alleles
    are partial junction sequences that do not match the full SV alleles in the
    truth VCF. Position-based matching allows monitoring of known SV loci.

    Args:
        sample_name: For logging.
        r1: R1 FASTQ path.
        r2: R2 FASTQ path.
        targets: Targets FASTA path.
        truth_vcf: Truth VCF for --target-variants.
        out_dir: Output directory.
        vtype: Variant type (used to decide whether to add position tolerance).

    Returns:
        True on success.
    """
    tmp = out_dir / "tmp_ti"
    tmp.mkdir(parents=True, exist_ok=True)

    cmd = [
        str(KAM), "run",
        "--r1", str(r1),
        "--r2", str(r2),
        "--targets", str(targets),
        "--target-variants", str(truth_vcf),
        "--output-dir", str(tmp),
        "--output-format", "vcf,tsv",
    ]
    if vtype in SV_TYPES:
        cmd += ["--ti-position-tolerance", "10"]

    result = subprocess.run(
        cmd,
        cwd=str(REPO),
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        print(
            f"[KAM TI FAIL] {sample_name}: {result.stderr[-600:]}",
            file=sys.stderr,
            flush=True,
        )
        shutil.rmtree(tmp, ignore_errors=True)
        return False

    for src, dst_name in [
        (tmp / "variants.tsv", "tumour_informed.tsv"),
        (tmp / "variants.vcf", "tumour_informed.vcf"),
    ]:
        if src.exists():
            try:
                shutil.copy2(src, out_dir / dst_name)
            except OSError as exc:
                print(f"[WARN] copy failed {src}: {exc}", file=sys.stderr)

    shutil.rmtree(tmp, ignore_errors=True)
    return True


def write_params_json(out_dir: Path, sample: dict[str, Any]) -> None:
    """Write params.json to the sample output directory.

    Stores all varforge parameters so the feature extraction script can
    include simulation conditions as features without reading the YAML.

    Args:
        out_dir: Sample output directory.
        sample: Sample dict from the manifest.
    """
    params = {
        "name": sample["name"],
        "split": sample["split"],
        "vtype": sample["vtype"],
        "coverage": sample["coverage"],
        "family_size_mean": sample["family_size_mean"],
        "pcr_cycles": sample["pcr_cycles"],
        "fragment_mean": sample["fragment_mean"],
        "fragment_sd": sample.get("fragment_sd"),
        "mean_quality": sample.get("mean_quality"),
        "vaf": sample["vaf"],
        "purity": sample["purity"],
        "seed": sample["seed"],
    }
    with open(out_dir / "params.json", "w") as fh:
        json.dump(params, fh, indent=2)


def is_complete(sample_name: str, split: str) -> bool:
    """Return True if a sample already has all required outputs.

    Args:
        sample_name: Sample identifier.
        split: 'train' or 'test'.

    Returns:
        True if both discovery.tsv and tumour_informed.tsv exist.
    """
    out_dir = SIM_DIR / split / sample_name
    return (
        (out_dir / "discovery.tsv").exists()
        and (out_dir / "tumour_informed.tsv").exists()
    )


# ─── Per-sample worker ────────────────────────────────────────────────────────

def process_sample(sample: dict[str, Any], dry_run: bool = False) -> tuple[str, str]:
    """Run the full pipeline for a single sample.

    Steps: varforge → kam discovery → kam tumour-informed → params.json.

    Args:
        sample: Sample dict from the manifest.
        dry_run: If True, skip actual execution.

    Returns:
        Tuple (sample_name, status) where status is 'done', 'skipped', or
        one of several failure strings.
    """
    name = sample["name"]
    split = sample["split"]
    vtype = sample["vtype"]
    out_dir = SIM_DIR / split / name

    if is_complete(name, split):
        return name, "skipped"

    if dry_run:
        return name, "dry_run"

    out_dir.mkdir(parents=True, exist_ok=True)

    # Determine config path
    config_dir = CONFIGS_TRAIN if split == "train" else CONFIGS_TEST
    config_path = config_dir / f"{name}.yaml"
    if not config_path.exists():
        print(f"[ERROR] Config not found: {config_path}", file=sys.stderr)
        return name, "no_config"

    # Step 1: varforge simulation
    if not any(out_dir.glob("*_R1.fastq.gz")):
        if not run_varforge(config_path, name):
            return name, "varforge_failed"

    # Step 2: locate FASTQs
    r1, r2 = find_fastqs(out_dir)
    if r1 is None:
        print(f"[ERROR] No FASTQs after simulation: {out_dir}", file=sys.stderr)
        return name, "no_fastq"

    # Step 3: locate truth VCF
    truth_vcf = find_truth_vcf(out_dir)
    if truth_vcf is None:
        # Fall back to the pre-generated truth VCF
        truth_vcf = ML_DIR / "truth_vcfs" / f"{name}.vcf"
        if not truth_vcf.exists():
            print(f"[ERROR] No truth VCF for {name}", file=sys.stderr)
            return name, "no_truth_vcf"

    targets = targets_for(vtype)

    # Step 4: kam discovery
    if not (out_dir / "discovery.tsv").exists():
        if not run_kam_discovery(name, r1, r2, targets, out_dir):
            return name, "disc_failed"

    # Step 5: kam tumour-informed
    if not (out_dir / "tumour_informed.tsv").exists():
        if not run_kam_tumour_informed(name, r1, r2, targets, truth_vcf, out_dir, vtype):
            return name, "ti_failed"

    # Step 6: write params
    write_params_json(out_dir, sample)

    print(f"[DONE] {name}", flush=True)
    return name, "done"


# ─── Nextcloud upload ─────────────────────────────────────────────────────────

def upload_batch(batch_samples: list[dict[str, Any]], batch_idx: int, split: str) -> bool:
    """Tar and upload a batch of sample directories to Nextcloud.

    Creates a tarball of all sample directories in the batch, uploads via
    WebDAV, then removes the local tarball.

    Args:
        batch_samples: List of sample dicts in this batch.
        batch_idx: Zero-based batch index.
        split: 'train' or 'test'.

    Returns:
        True on successful upload.
    """
    batch_name = f"batch_{batch_idx:03d}"
    tarball = Path(f"/tmp/ml_twist_{split}_{batch_name}.tar.gz")
    split_dir = SIM_DIR / split

    # Collect sample directory names that actually exist
    dirs_to_tar = [
        s["name"]
        for s in batch_samples
        if (split_dir / s["name"]).is_dir()
    ]

    if not dirs_to_tar:
        print(f"[UPLOAD] No directories to upload for {batch_name}", flush=True)
        return True

    # Build tarball
    tar_cmd = ["tar", "-czf", str(tarball), "-C", str(split_dir)] + dirs_to_tar
    tar_result = subprocess.run(tar_cmd, capture_output=True, text=True)
    if tar_result.returncode != 0:
        print(
            f"[UPLOAD FAIL] tar failed for {batch_name}: {tar_result.stderr[-400:]}",
            file=sys.stderr,
        )
        return False

    # Upload via WebDAV
    upload_url = (
        f"{NEXTCLOUD_URL}/benchmarking/ml-twist-duplex/{split}/{batch_name}.tar.gz"
    )
    curl_result = subprocess.run(
        [
            "curl", "--silent", "--show-error",
            "-T", str(tarball),
            upload_url,
            "-u", f"{NEXTCLOUD_TOKEN}:",
        ],
        capture_output=True,
        text=True,
    )

    try:
        tarball.unlink()
    except OSError:
        pass

    if curl_result.returncode != 0:
        print(
            f"[UPLOAD FAIL] curl failed for {batch_name}: {curl_result.stderr[-400:]}",
            file=sys.stderr,
        )
        return False

    print(f"[UPLOAD] {split}/{batch_name} uploaded ({len(dirs_to_tar)} dirs)", flush=True)
    return True


# ─── Batch processing ────────────────────────────────────────────────────────

def process_batch(
    batch: list[dict[str, Any]],
    batch_idx: int,
    checkpoint: dict[str, str],
    workers: int,
    dry_run: bool,
    skip_upload: bool,
) -> dict[str, str]:
    """Process a batch of samples in parallel and upload results.

    Args:
        batch: List of sample dicts.
        batch_idx: Zero-based batch index.
        checkpoint: Mutable checkpoint dict (updated in place).
        workers: Number of parallel workers.
        dry_run: If True, skip execution and upload.
        skip_upload: If True, skip the Nextcloud upload step.

    Returns:
        Updated checkpoint dict.
    """
    split = batch[0]["split"]
    to_run = [s for s in batch if checkpoint.get(s["name"]) != "done"]

    if not to_run:
        print(f"[BATCH {batch_idx}] All {len(batch)} samples already done — skipping.", flush=True)
        return checkpoint

    print(f"\n[BATCH {batch_idx}] Processing {len(to_run)} samples (split={split})...", flush=True)

    done = skipped = failed = 0
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(process_sample, s, dry_run): s for s in to_run}
        for fut in as_completed(futures):
            sample = futures[fut]
            try:
                name, status = fut.result()
            except Exception as exc:
                name = sample["name"]
                status = f"exception:{exc}"
                print(f"[ERROR] {name}: {exc}", file=sys.stderr, flush=True)

            checkpoint[name] = status
            if status == "done":
                done += 1
            elif status in ("skipped", "dry_run"):
                skipped += 1
            else:
                failed += 1

    save_checkpoint(checkpoint)
    print(
        f"[BATCH {batch_idx}] Done={done} Skipped={skipped} Failed={failed}",
        flush=True,
    )

    if not dry_run and not skip_upload:
        upload_batch(batch, batch_idx, split)

    return checkpoint


# ─── CLI ──────────────────────────────────────────────────────────────────────

def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Run the full Twist duplex ML pipeline (varforge + kam).",
    )
    parser.add_argument(
        "--split",
        choices=["train", "test", "all"],
        default="all",
        help="Which split to process (default: all).",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=DEFAULT_WORKERS,
        help=f"Number of parallel workers (default: {DEFAULT_WORKERS}).",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=DEFAULT_BATCH_SIZE,
        help=f"Samples per batch (default: {DEFAULT_BATCH_SIZE}).",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print what would run without executing anything.",
    )
    parser.add_argument(
        "--skip-upload",
        action="store_true",
        help="Skip Nextcloud upload after each batch.",
    )
    return parser.parse_args()


def main() -> None:
    """Entry point."""
    args = parse_args()

    if not KAM.exists():
        print(f"[ERROR] kam binary not found: {KAM}", file=sys.stderr)
        print("Build with: cargo build --release", file=sys.stderr)
        sys.exit(1)

    samples = load_manifest()

    # Filter by split
    if args.split != "all":
        samples = [s for s in samples if s["split"] == args.split]

    print(f"Total samples to process: {len(samples)}", flush=True)

    checkpoint = load_checkpoint()
    already_done = sum(1 for s in samples if checkpoint.get(s["name"]) == "done")
    print(f"Already completed: {already_done}", flush=True)

    # Split into batches
    batches = [
        samples[i : i + args.batch_size]
        for i in range(0, len(samples), args.batch_size)
    ]
    print(f"Batches: {len(batches)} x ~{args.batch_size} samples", flush=True)

    t_start = time.time()

    for batch_idx, batch in enumerate(batches):
        checkpoint = process_batch(
            batch,
            batch_idx,
            checkpoint,
            workers=args.workers,
            dry_run=args.dry_run,
            skip_upload=args.skip_upload,
        )

    elapsed = time.time() - t_start
    total_done = sum(1 for v in checkpoint.values() if v == "done")
    total_failed = sum(1 for v in checkpoint.values() if v not in ("done", "skipped", "dry_run"))
    print(f"\n=== COMPLETE ===")
    print(f"Done:   {total_done}")
    print(f"Failed: {total_failed}")
    print(f"Elapsed: {elapsed / 60:.1f} min")


if __name__ == "__main__":
    main()
