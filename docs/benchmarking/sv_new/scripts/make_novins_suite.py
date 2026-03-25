#!/usr/bin/env python3
"""Generate varforge configs and truth VCFs for NovelInsertion benchmarking.

Produces 25 VAF levels × 2 replicates = 50 configs and 50 truth VCFs.
The inserted sequence is 80 bp of non-repetitive synthetic DNA that does not
appear (as any 20-mer) in the 2000 bp reference or its reverse complement.
It is also not a tandem duplication of any reference window.

The canonical insertion sequence is also written to:
    docs/benchmarking/sv_new/data/novins_insert_seq.txt

Insertion site: chr1:400 (same position as the existing INS suite).
SVTYPE=INS; the NovelInsertion classifier distinguishes this from tandem
duplications by verifying the inserted sequence is absent from the reference.

Run from the repo root:
    python3 docs/benchmarking/sv_new/scripts/make_novins_suite.py
"""

from pathlib import Path

VAF_LEVELS = [
    0.0005, 0.0010, 0.0015, 0.0020, 0.0025, 0.0030, 0.0035, 0.0040, 0.0050,
    0.0060, 0.0075, 0.0100, 0.0125, 0.0150, 0.0175, 0.0200, 0.0250, 0.0300,
    0.0350, 0.0400, 0.0500, 0.0600, 0.0700, 0.0800, 0.1000,
]

ROOT = Path("docs/benchmarking/sv_new")
DATA = ROOT / "data"
CFGS = ROOT / "configs"

# Canonical novel insertion sequence: 80 bp, generated deterministically from
# seed 202503251 (random.seed). Verified: no 20-mer appears in the 2000 bp
# synthetic reference or its reverse complement.
NOVEL_INSERT_SEQ = (
    "CCGGCGTCACCCGATGATCAGGTGCCTAGAATCAAGCGTAAGGGCTCGTCGTCACATGTGATAATCTTCCGTGTCCAGCC"
)

# Reference base at chr1:400 (1-based). Position 400 in the VCF is 1-based
# so the REF allele is a single base. The ALT is REF + inserted sequence.
# From ref.fa: position 399 (0-based) in the concatenated sequence.
REF_BASE_AT_400 = "T"

NOVINS_HEADERS = [
    "##fileformat=VCFv4.2",
    "##contig=<ID=chr1,length=2000>",
    '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV type">',
    '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="SV length (positive=insertion)">',
    '##INFO=<ID=SEQ,Number=1,Type=String,Description="Inserted sequence">',
    '##INFO=<ID=VAF,Number=1,Type=Float,Description="Variant allele frequency">',
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
]


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


def vaf_tag(vaf: float) -> str:
    """Return zero-padded 4-digit VAF tag (e.g. 0.01 → '0100')."""
    return f"{round(vaf * 10000):04d}"


def write_vcf(path: Path, vaf: float) -> None:
    """Write a truth VCF for one NovelInsertion sample at the given VAF."""
    alt = REF_BASE_AT_400 + NOVEL_INSERT_SEQ
    svlen = len(NOVEL_INSERT_SEQ)
    row = (
        f"chr1\t400\t.\t{REF_BASE_AT_400}\t{alt}\t.\tPASS"
        f"\tSVTYPE=INS;SVLEN={svlen};SEQ={NOVEL_INSERT_SEQ};VAF={vaf}"
    )
    with open(path, "w") as f:
        for h in NOVINS_HEADERS:
            f.write(h + "\n")
        f.write(row + "\n")


def write_config(
    path: Path,
    sample_name: str,
    vaf: float,
    purity: float,
    vcf_rel: str,
    out_dir: str,
    seed: int,
) -> None:
    """Write a varforge YAML config for one NovelInsertion sample."""
    with open(path, "w") as f:
        f.write(CONFIG_TEMPLATE.format(
            sample_name=sample_name,
            vaf_pct=vaf * 100,
            purity=purity,
            vcf_path=vcf_rel,
            out_dir=out_dir,
            seed=seed,
        ))


def write_insert_seq_file() -> None:
    """Write the canonical inserted sequence to a reference file."""
    out = DATA / "novins_insert_seq.txt"
    with open(out, "w") as f:
        f.write("# Canonical novel insertion sequence for NovelInsertion benchmarking.\n")
        f.write("# Length: 80 bp. Seed: 202503251 (Python random).\n")
        f.write("# Verified: no 20-mer appears in the 2000 bp synthetic reference\n")
        f.write("# (docs/benchmarking/sv/data/ref.fa) or its reverse complement.\n")
        f.write(NOVEL_INSERT_SEQ + "\n")
    print(f"Wrote canonical insert sequence to {out}")


def main() -> None:
    DATA.mkdir(parents=True, exist_ok=True)
    CFGS.mkdir(parents=True, exist_ok=True)

    write_insert_seq_file()

    count = 0
    for vaf in VAF_LEVELS:
        t = vaf_tag(vaf)
        purity = round(vaf * 2, 6)
        tag_int = round(vaf * 10000)

        for rep_idx, rep in enumerate(["a", "b"]):
            rep_up = rep.upper()
            seed_offset = rep_idx * 1000

            sample_name = f"NOVINS_VAF{t}_{rep_up}"
            vcf_path = DATA / f"truth_novins_vaf{t}_{rep}.vcf"
            cfg_path = CFGS / f"novins_vaf{t}_{rep}.yaml"
            vcf_rel = f"docs/benchmarking/sv_new/data/truth_novins_vaf{t}_{rep}.vcf"
            out_dir = f"docs/benchmarking/sv_new/results/sim_novins_vaf{t}_{rep}"
            # Seed range 110000+ avoids collision with all existing suites
            # (sv: 50000+, ins: 70000+, invdel: 90000+, sv_new invdel: 90000+).
            seed = 110000 + tag_int + seed_offset

            write_vcf(vcf_path, vaf)
            write_config(cfg_path, sample_name, vaf, purity, vcf_rel, out_dir, seed)
            count += 1

    print(f"Generated {count} NovelInsertion configs and truth VCFs in {ROOT}/")


if __name__ == "__main__":
    main()
