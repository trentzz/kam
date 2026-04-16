#!/usr/bin/env python3
"""Generate varforge configs for the SV size sweep experiment.

Creates configs for DEL, DUP, and INV at multiple SV sizes (20-500 bp).
INV configs are only generated for sizes >= 50 bp since smaller inversions
are classified as indels.

Each config follows the 03-sv-extended format with inline: false UMIs.

Output goes to scripts/benchmarking/sv_size_sweep_configs/.
"""

import os
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
CONFIGS_DIR = SCRIPT_DIR / "sv_size_sweep_configs"

# SV sizes to sweep.
SIZES = [20, 50, 100, 200, 500]

# SV types and their size constraints.
SV_TYPES = {
    "del": {"min_size": 20},
    "dup": {"min_size": 20},
    "inv": {"min_size": 50},  # inversions < 50 bp are indels
}

# Two replicates per config.
REPLICATES = ["a", "b"]

# Seed ranges: 200000+ for size sweep, offset by type and size.
SEED_BASE = {
    "del": 200000,
    "dup": 210000,
    "inv": 220000,
}

TEMPLATE = """\
# varforge config: {name} — {size} bp {sv_type_upper} at VAF 1.000%
reference: docs/benchmarking/sv/data/ref.fa

output:
  directory: docs/benchmarking/sv_new/results/sim_{sv_type}_{size}bp_vaf0100_{rep}
  fastq: true
  bam: false
  truth_vcf: true
  manifest: true

sample:
  name: {name}
  read_length: 150
  coverage: 5000.0
  platform: illumina

fragment:
  model: normal
  mean: 167.0
  sd: 30.0

quality:
  mean_quality: 37
  tail_decay: 0.002

tumour:
  purity: 0.020000
  ploidy: 2

mutations:
  vcf: docs/benchmarking/sv_new/data/truth_{sv_type}_{size}bp_vaf0100_{rep}.vcf

umi:
  length: 5
  duplex: true
  pcr_cycles: 8
  family_size_mean: 4.0
  family_size_sd: 1.5
  inline: false

chromosomes:
  - chr1

seed: {seed}
"""


def main() -> None:
    CONFIGS_DIR.mkdir(parents=True, exist_ok=True)

    count = 0
    for sv_type, constraints in SV_TYPES.items():
        min_size = constraints["min_size"]
        for size in SIZES:
            if size < min_size:
                continue
            for i, rep in enumerate(REPLICATES):
                name = f"{sv_type.upper()}_{size}BP_VAF0100_{rep.upper()}"
                seed = SEED_BASE[sv_type] + size * 10 + i

                content = TEMPLATE.format(
                    name=name,
                    size=size,
                    sv_type=sv_type,
                    sv_type_upper=sv_type.upper(),
                    rep=rep,
                    seed=seed,
                )

                filename = f"{sv_type}_{size}bp_vaf0100_{rep}.yaml"
                cfg_path = CONFIGS_DIR / filename
                with open(cfg_path, "w") as f:
                    f.write(content)
                count += 1

    print(f"Generated {count} size sweep configs in {CONFIGS_DIR}")


if __name__ == "__main__":
    main()
