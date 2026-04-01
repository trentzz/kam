#!/usr/bin/env python3
"""Generate varforge configs for the Twist duplex ML dataset.

Produces 11,000 unique samples (10,000 train + 1,000 test) by selecting from
existing truth VCF pools. Each sample gets a unique YAML config with randomised
simulation parameters.

VCF pools used:
  docs/benchmarking/snvindel/data/truth_snvs_vaf*.vcf      (54 files)
  docs/benchmarking/snvindel/data/truth_indels_vaf*.vcf    (54 files)
  docs/benchmarking/sv/data/truth_svs_vaf*.vcf             (55 files)
  docs/benchmarking/sv/data/truth_ins_vaf*.vcf             (54 files)
  docs/benchmarking/sv/data/truth_invdel_vaf*.vcf          (54 files)

Sample counts:
  snv:       3,500 train + 350 test
  indel:     3,000 train + 300 test
  sv_dupinv: 1,500 train + 150 test
  ins:       1,500 train + 150 test
  invdel:      500 train +  50 test

Usage:
    python3 scripts/ml/generate_twist_duplex_configs.py [--dry-run]

Run from the repository root.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

import numpy as np

REPO = Path(__file__).resolve().parent.parent.parent

# ─── Output paths ─────────────────────────────────────────────────────────────

ML_DIR = REPO / "docs" / "benchmarking" / "ml-twist-duplex"
CONFIGS_DIR = ML_DIR / "configs"
SIMS_DIR = ML_DIR / "simulations"
MANIFEST_PATH = ML_DIR / "manifest.json"

# ─── Reference files (relative to REPO, as used in YAML) ─────────────────────

SNVINDEL_REF = "docs/benchmarking/snvindel/data/ref.fa"
SV_REF = "docs/benchmarking/sv/data/ref.fa"

# ─── VCF pool directories ─────────────────────────────────────────────────────

SNVINDEL_DATA = REPO / "docs" / "benchmarking" / "snvindel" / "data"
SV_DATA = REPO / "docs" / "benchmarking" / "sv" / "data"

# ─── Targets FASTA (relative to REPO, stored in manifest) ────────────────────

SNVINDEL_TARGETS = "docs/benchmarking/snvindel/data/snvindel_targets.fa"
SV_SUITE_TARGETS = "docs/benchmarking/sv/data/sv_suite_targets.fa"
INS_TARGETS = "docs/benchmarking/sv/data/ins_targets.fa"
INVDEL_TARGETS = "docs/benchmarking/sv/data/invdel_targets.fa"

# ─── Dataset spec ─────────────────────────────────────────────────────────────

# (vtype, n_train, n_test, vcf_glob, reference, targets)
VARIANT_TYPES: list[tuple[str, int, int, str, str, str]] = [
    ("snv",       3500, 350, "truth_snvs_vaf*.vcf",    SNVINDEL_REF, SNVINDEL_TARGETS),
    ("indel",     3000, 300, "truth_indels_vaf*.vcf",  SNVINDEL_REF, SNVINDEL_TARGETS),
    ("sv_dupinv", 1500, 150, "truth_svs_vaf*.vcf",     SV_REF,       SV_SUITE_TARGETS),
    ("ins",       1500, 150, "truth_ins_vaf*.vcf",     SV_REF,       INS_TARGETS),
    ("invdel",     500,  50, "truth_invdel_vaf*.vcf",  SV_REF,       INVDEL_TARGETS),
]

# VCF pool directories: snvindel types use SNVINDEL_DATA, SV types use SV_DATA
VCF_POOL_DIR: dict[str, Path] = {
    "snv":       SNVINDEL_DATA,
    "indel":     SNVINDEL_DATA,
    "sv_dupinv": SV_DATA,
    "ins":       SV_DATA,
    "invdel":    SV_DATA,
}


# ─── RNG seed formula ────────────────────────────────────────────────────────

def sample_seed(global_idx: int) -> int:
    """Return a deterministic seed for the given global sample index.

    Args:
        global_idx: Zero-based index across all 11,000 samples.

    Returns:
        Integer seed for varforge.
    """
    return global_idx * 9973 + 42


# ─── Parameter sampling ───────────────────────────────────────────────────────

def sample_params(rng: np.random.Generator) -> dict[str, Any]:
    """Sample randomised simulation parameters for one sample.

    All parameters are drawn independently per sample to maximise
    diversity in the training set.

    Args:
        rng: NumPy random generator.

    Returns:
        Dict with keys: vaf, purity, coverage, family_size_mean,
        family_size_sd, pcr_cycles, fragment_mean, fragment_sd,
        mean_quality.
    """
    # VAF: log-uniform 0.01% – 5%
    log_vaf = rng.uniform(np.log(0.0001), np.log(0.05))
    vaf = float(np.exp(log_vaf))
    purity = vaf * 2.0  # diploid: purity = 2 * VAF

    # Coverage: log-uniform 1,000 – 15,000
    log_cov = rng.uniform(np.log(1000.0), np.log(15000.0))
    coverage = float(np.exp(log_cov))

    family_size_mean = float(rng.uniform(2.0, 8.0))
    pcr_cycles = int(rng.integers(5, 15))  # 5 inclusive, 14 inclusive
    fragment_mean = float(rng.uniform(130.0, 280.0))
    fragment_sd = float(rng.uniform(20.0, 50.0))
    mean_quality = int(rng.integers(30, 41))  # 30–40 inclusive

    return {
        "vaf": vaf,
        "purity": purity,
        "coverage": coverage,
        "family_size_mean": family_size_mean,
        "family_size_sd": 1.0,
        "pcr_cycles": pcr_cycles,
        "fragment_mean": fragment_mean,
        "fragment_sd": fragment_sd,
        "mean_quality": mean_quality,
    }


# ─── YAML generation ─────────────────────────────────────────────────────────

def build_yaml(
    sample_name: str,
    split: str,
    vtype: str,
    reference: str,
    truth_vcf_rel: str,
    params: dict[str, Any],
    seed: int,
) -> str:
    """Render a varforge config YAML string.

    Args:
        sample_name: Unique sample identifier, e.g. 'snv_00001'.
        split: 'train' or 'test'.
        vtype: Variant type string.
        reference: Relative path to the reference FASTA.
        truth_vcf_rel: Relative path from REPO root to the chosen truth VCF.
        params: Dict from sample_params().
        seed: Integer seed.

    Returns:
        YAML string ready to write to disk.
    """
    vaf_pct = params["vaf"] * 100.0
    sim_dir = (
        f"docs/benchmarking/ml-twist-duplex/simulations/{split}/{sample_name}"
    )

    lines = [
        f"# varforge config: {sample_name}  — VAF {vaf_pct:.4f}%  type={vtype}",
        f"reference: {reference}",
        "",
        "output:",
        f"  directory: {sim_dir}",
        "  fastq: true",
        "  bam: false",
        "  truth_vcf: true",
        "  manifest: true",
        "",
        "sample:",
        f"  name: {sample_name}",
        "  read_length: 150",
        f"  coverage: {params['coverage']:.1f}",
        "  platform: illumina",
        "",
        "fragment:",
        "  model: normal",
        f"  mean: {params['fragment_mean']:.1f}",
        f"  sd: {params['fragment_sd']:.1f}",
        "",
        "quality:",
        f"  mean_quality: {params['mean_quality']}",
        "  tail_decay: 0.002",
        "",
        "tumour:",
        f"  purity: {params['purity']:.6f}",
        "  ploidy: 2",
        "",
        "mutations:",
        f"  vcf: {truth_vcf_rel}",
        "",
        "umi:",
        "  length: 5",
        "  duplex: true",
        f"  pcr_cycles: {params['pcr_cycles']}",
        f"  family_size_mean: {params['family_size_mean']:.2f}",
        f"  family_size_sd: {params['family_size_sd']:.1f}",
        "  inline: true",
        "",
        "chromosomes:",
        "  - chr1",
        "",
        f"seed: {seed}",
        "",
    ]
    return "\n".join(lines)


# ─── Sample name helpers ──────────────────────────────────────────────────────

def make_name(vtype: str, local_idx: int, split: str) -> str:
    """Build the sample name for a given type, local index, and split.

    Train samples: {type}_{idx:05d}  e.g. snv_00001
    Test  samples: {type}_t{idx:04d} e.g. snv_t0001

    Local index is 1-based.

    Args:
        vtype: Variant type string.
        local_idx: 1-based index within this type+split.
        split: 'train' or 'test'.

    Returns:
        Sample name string.
    """
    if split == "train":
        return f"{vtype}_{local_idx:05d}"
    return f"{vtype}_t{local_idx:04d}"


# ─── VCF pool loading ─────────────────────────────────────────────────────────

def load_vcf_pool(vtype: str, glob_pattern: str) -> list[Path]:
    """Return all VCF paths matching the glob in the pool directory.

    Args:
        vtype: Variant type (used to select the pool directory).
        glob_pattern: Glob pattern, e.g. 'truth_snvs_vaf*.vcf'.

    Returns:
        Sorted list of matching paths.

    Raises:
        SystemExit: If no VCF files are found.
    """
    pool_dir = VCF_POOL_DIR[vtype]
    pool = sorted(pool_dir.glob(glob_pattern))
    if not pool:
        print(
            f"[ERROR] No VCF files matching '{glob_pattern}' in {pool_dir}",
            file=sys.stderr,
        )
        sys.exit(1)
    return pool


# ─── Main ─────────────────────────────────────────────────────────────────────

def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Generate varforge configs for the Twist duplex ML dataset.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print summary without writing any files.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=12345,
        help="Global RNG seed for parameter sampling (default: 12345).",
    )
    return parser.parse_args()


def main() -> None:
    """Generate all configs and write manifest.json."""
    args = parse_args()

    rng = np.random.default_rng(args.seed)

    # Load VCF pools once per variant type
    vcf_pools: dict[str, list[Path]] = {}
    for vtype, _n_train, _n_test, glob_pat, _ref, _tgt in VARIANT_TYPES:
        vcf_pools[vtype] = load_vcf_pool(vtype, glob_pat)
        print(
            f"  Pool [{vtype:10s}]: {len(vcf_pools[vtype])} VCFs",
            flush=True,
        )

    if not args.dry_run:
        for split in ("train", "test"):
            (CONFIGS_DIR / split).mkdir(parents=True, exist_ok=True)
            (SIMS_DIR / split).mkdir(parents=True, exist_ok=True)
        ML_DIR.mkdir(parents=True, exist_ok=True)

    all_samples: list[dict[str, Any]] = []
    global_idx = 0
    total_written = 0

    for vtype, n_train, n_test, _glob, reference, targets in VARIANT_TYPES:
        pool = vcf_pools[vtype]

        for split, n_samples in (("train", n_train), ("test", n_test)):
            # Sample VCF indices with replacement; pool is small, samples large
            vcf_indices = rng.integers(0, len(pool), size=n_samples)

            for local_idx in range(1, n_samples + 1):
                seed = sample_seed(global_idx)
                params = sample_params(rng)

                name = make_name(vtype, local_idx, split)
                chosen_vcf: Path = pool[int(vcf_indices[local_idx - 1])]
                # Store relative to REPO for portability
                truth_vcf_rel = str(chosen_vcf.relative_to(REPO))

                yaml_str = build_yaml(
                    sample_name=name,
                    split=split,
                    vtype=vtype,
                    reference=reference,
                    truth_vcf_rel=truth_vcf_rel,
                    params=params,
                    seed=seed,
                )

                config_path = CONFIGS_DIR / split / f"{name}.yaml"

                if not args.dry_run:
                    config_path.write_text(yaml_str)
                    total_written += 1

                all_samples.append({
                    "name": name,
                    "split": split,
                    "vtype": vtype,
                    "seed": seed,
                    "vaf": round(params["vaf"], 8),
                    "purity": round(params["purity"], 8),
                    "coverage": round(params["coverage"], 2),
                    "family_size_mean": round(params["family_size_mean"], 4),
                    "family_size_sd": params["family_size_sd"],
                    "pcr_cycles": params["pcr_cycles"],
                    "fragment_mean": round(params["fragment_mean"], 2),
                    "fragment_sd": round(params["fragment_sd"], 2),
                    "mean_quality": params["mean_quality"],
                    "truth_vcf": truth_vcf_rel,
                    "config": str(
                        (CONFIGS_DIR / split / f"{name}.yaml").relative_to(REPO)
                    ),
                    "reference": reference,
                    "targets": targets,
                    "sim_dir": (
                        f"docs/benchmarking/ml-twist-duplex/simulations/"
                        f"{split}/{name}"
                    ),
                })

                global_idx += 1

    # Write manifest
    manifest: dict[str, Any] = {
        "version": 1,
        "n_train": sum(1 for s in all_samples if s["split"] == "train"),
        "n_test": sum(1 for s in all_samples if s["split"] == "test"),
        "samples": all_samples,
    }

    if not args.dry_run:
        with open(MANIFEST_PATH, "w") as fh:
            json.dump(manifest, fh, indent=2)
        print(f"\nManifest: {MANIFEST_PATH}", flush=True)
        print(f"Configs written: {total_written}", flush=True)
    else:
        print(f"\n[DRY RUN] Would write {len(all_samples)} configs.", flush=True)

    # Summary
    print(f"\nTotal samples: {len(all_samples)}", flush=True)
    for vtype, n_train, n_test, _glob, _ref, _tgt in VARIANT_TYPES:
        print(
            f"  {vtype:10s}: {n_train} train + {n_test} test "
            f"({n_train + n_test} total)",
            flush=True,
        )


if __name__ == "__main__":
    main()
