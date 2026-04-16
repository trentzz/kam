#!/usr/bin/env python3
"""Generate varforge configs for ultra-low VAF SV detection limits.

Creates configs at 0.01%, 0.02%, 0.03%, 0.04% VAF for each SV type
(DEL, DUP, INV, INS, InvDel, NovIns) at 5000x coverage.

Each config follows the 03-sv-extended format with inline: false UMIs.

Output goes to docs/benchmarking/03-sv-extended/configs/ultra_low_vaf/.
"""

import os
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[1]
CONFIGS_DIR = REPO_ROOT / "docs" / "benchmarking" / "03-sv-extended" / "configs" / "ultra_low_vaf"

# Ultra-low VAF levels (as fractions).
VAF_LEVELS = [
    (0.0001, "0001", "0.010"),
    (0.0002, "0002", "0.020"),
    (0.0003, "0003", "0.030"),
    (0.0004, "0004", "0.040"),
]

# SV types to generate configs for.
SV_TYPES = ["del", "dup", "inv", "ins", "invdel", "novins"]

# Two replicates.
REPLICATES = ["a", "b"]

# Seed base per type (300000+ range to avoid collisions with other suites).
SEED_BASE = {
    "del": 300000,
    "dup": 310000,
    "inv": 320000,
    "ins": 330000,
    "invdel": 340000,
    "novins": 350000,
}

TEMPLATE = """\
# varforge config: {name} — {sv_type_upper} at VAF {vaf_pct}%
reference: docs/benchmarking/sv/data/ref.fa

output:
  directory: docs/benchmarking/sv_new/results/sim_{sv_type}_vaf{vaf_tag}_{rep}
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
  purity: {purity:.6f}
  ploidy: 2

mutations:
  vcf: docs/benchmarking/sv_new/data/truth_{sv_type}_vaf{vaf_tag}_{rep}.vcf

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
    for sv_type in SV_TYPES:
        for vaf, vaf_tag, vaf_pct in VAF_LEVELS:
            purity = vaf * 2  # diploid het: effective VAF = purity / 2
            for i, rep in enumerate(REPLICATES):
                name = f"{sv_type.upper()}_VAF{vaf_tag}_{rep.upper()}"
                seed = SEED_BASE[sv_type] + int(vaf_tag) * 10 + i

                content = TEMPLATE.format(
                    name=name,
                    sv_type=sv_type,
                    sv_type_upper=sv_type.upper(),
                    vaf_tag=vaf_tag,
                    vaf_pct=vaf_pct,
                    purity=purity,
                    rep=rep,
                    seed=seed,
                )

                filename = f"{sv_type}_vaf{vaf_tag}_{rep}.yaml"
                cfg_path = CONFIGS_DIR / filename
                with open(cfg_path, "w") as f:
                    f.write(content)
                count += 1

    print(f"Generated {count} ultra-low VAF configs in {CONFIGS_DIR}")


if __name__ == "__main__":
    main()
