#!/usr/bin/env python3
"""Generate varforge configs and truth VCFs for the SV benchmark suite.

Produces 25 VAF levels × 2 replicates × 3 types (DUP+INV, INS, INVDEL)
= 150 configs + 150 truth VCFs.

Note: DEL variants (≥20bp) are excluded from the DUP+INV suite due to a
varforge engine panic (range bounds error in src/core/engine.rs:393 when
applying large deletions). DUP and INV are unaffected. Run from the repo root.
"""
import re
from pathlib import Path

VAF_LEVELS = [
    0.0005, 0.0010, 0.0015, 0.0020, 0.0025, 0.0030, 0.0035, 0.0040, 0.0050,
    0.0060, 0.0075, 0.0100, 0.0125, 0.0150, 0.0175, 0.0200, 0.0250, 0.0300,
    0.0350, 0.0400, 0.0500, 0.0600, 0.0700, 0.0800, 0.1000,
]

ROOT  = Path("docs/benchmarking/sv")
DATA  = ROOT / "data"
CFGS  = ROOT / "configs"

# Read variant rows from an existing truth VCF, stripping the VAF field.
def read_vcf(path: Path):
    headers, rows = [], []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("#"):
                headers.append(line)
            else:
                rows.append(line)
    return headers, rows


def set_vaf(row: str, vaf: float) -> str:
    return re.sub(r"VAF=[\d.]+", f"VAF={vaf}", row)


def write_vcf(path: Path, headers: list, rows: list, vaf: float) -> None:
    with open(path, "w") as f:
        for h in headers:
            f.write(h + "\n")
        for row in rows:
            f.write(set_vaf(row, vaf) + "\n")


CONFIG_TEMPLATE = """\
# varforge config: {sample_name} — VAF {vaf_pct:.3f}%
reference: docs/benchmarking/sv/data/ref.fa

output:
  directory: {out_dir}
  fastq: true
  bam: false
  truth_vcf: true
  manifest: true

sample:
  name: {sample_name}
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
  vcf: {vcf_path}

umi:
  length: 5
  duplex: true
  pcr_cycles: 8
  family_size_mean: 4.0
  family_size_sd: 1.5
  inline: true

chromosomes:
  - chr1

seed: {seed}
"""


def write_config(path: Path, sample_name: str, vaf: float, purity: float,
                 vcf_rel: str, out_dir: str, seed: int) -> None:
    with open(path, "w") as f:
        f.write(CONFIG_TEMPLATE.format(
            sample_name=sample_name,
            vaf_pct=vaf * 100,
            purity=purity,
            vcf_path=vcf_rel,
            out_dir=out_dir,
            seed=seed,
        ))


def vaf_tag(vaf: float) -> str:
    return f"{round(vaf * 10000):04d}"


# Load template rows from existing 1% truth VCFs.
sv_hdr,      sv_rows      = read_vcf(DATA / "truth_svs_vaf010.vcf")
ins_hdr,     ins_rows     = read_vcf(DATA / "truth_ins_vaf010.vcf")
invdel_hdr,  invdel_rows  = read_vcf(DATA / "truth_invdel_vaf010.vcf")

counts = {"sv": 0, "ins": 0, "invdel": 0}

for vaf in VAF_LEVELS:
    t = vaf_tag(vaf)
    purity  = round(vaf * 2, 6)
    tag_int = round(vaf * 10000)

    for rep_idx, rep in enumerate(["a", "b"]):
        rep_up = rep.upper()
        seed_offset = rep_idx * 1000

        # --- DUP+INV (DEL excluded: varforge panic for deletions ≥20bp) ---
        sv_vcf  = DATA / f"truth_svs_vaf{t}_{rep}.vcf"
        sv_cfg  = CFGS / f"sv_vaf{t}_{rep}.yaml"
        write_vcf(sv_vcf, sv_hdr, sv_rows, vaf)
        write_config(sv_cfg,
            sample_name=f"SV_VAF{t}_{rep_up}",
            vaf=vaf, purity=purity,
            vcf_rel=f"docs/benchmarking/sv/data/truth_svs_vaf{t}_{rep}.vcf",
            out_dir=f"docs/benchmarking/sv/results/sim_sv_vaf{t}_{rep}",
            seed=50000 + tag_int + seed_offset,
        )
        counts["sv"] += 1

        # --- INS ---
        ins_vcf  = DATA / f"truth_ins_vaf{t}_{rep}.vcf"
        ins_cfg  = CFGS / f"ins_vaf{t}_{rep}.yaml"
        write_vcf(ins_vcf, ins_hdr, ins_rows, vaf)
        write_config(ins_cfg,
            sample_name=f"INS_VAF{t}_{rep_up}",
            vaf=vaf, purity=purity,
            vcf_rel=f"docs/benchmarking/sv/data/truth_ins_vaf{t}_{rep}.vcf",
            out_dir=f"docs/benchmarking/sv/results/sim_ins_vaf{t}_{rep}",
            seed=70000 + tag_int + seed_offset,
        )
        counts["ins"] += 1

        # --- INVDEL ---
        invdel_vcf = DATA / f"truth_invdel_vaf{t}_{rep}.vcf"
        invdel_cfg = CFGS / f"invdel_vaf{t}_{rep}.yaml"
        write_vcf(invdel_vcf, invdel_hdr, invdel_rows, vaf)
        write_config(invdel_cfg,
            sample_name=f"INVDEL_VAF{t}_{rep_up}",
            vaf=vaf, purity=purity,
            vcf_rel=f"docs/benchmarking/sv/data/truth_invdel_vaf{t}_{rep}.vcf",
            out_dir=f"docs/benchmarking/sv/results/sim_invdel_vaf{t}_{rep}",
            seed=90000 + tag_int + seed_offset,
        )
        counts["invdel"] += 1

total = sum(counts.values())
print(f"Generated {counts['sv']} SV + {counts['ins']} INS + {counts['invdel']} INVDEL configs ({total} total)")
