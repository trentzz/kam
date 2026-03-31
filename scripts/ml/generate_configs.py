"""Generate varforge configs for all ML sim dirs that don't have FASTQs.

Regenerates configs based on the patterns observed in existing manifest files.
"""

import re
from pathlib import Path

REPO = Path("/home/trent/code/kam")
CONFIGS_DIR = REPO / "docs/benchmarking/ml/configs"
RESULTS = REPO / "docs/benchmarking/ml/results"

# Seed base per variant type and replicate
# Pattern: type_base + rep_offset + vaf_int
TYPE_BASES = {
    "snv": 10000,
    "indel": 30000,
    "sv": 50000,
    "ins": 70000,
    "invdel": 90000,
}
REP_OFFSETS = {
    "a": 0,
    "b": 1000,
    "c": 2000,
    "d": 3000,
    "e": 4000,
}

# Parameter sweep values
COVERAGES = {
    "cov1000": 1000.0,
    "cov2000": 2000.0,
    "cov8000": 8000.0,
    "cov12000": 12000.0,
}
FAM_SIZES = {
    "fam2": 2.0,
    "fam6": 6.0,
    "fam8": 8.0,
}
PCR_CYCLES = {
    "pcr5": 5,
    "pcr10": 10,
    "pcr12": 12,
}
FRAG_MEANS = {
    "frag140": 140.0,
    "frag200": 200.0,
    "frag250": 250.0,
}

# Mutation VCF paths per type (always _a variant)
MUTATION_VCFS = {
    "snv": "docs/benchmarking/snvindel/data/truth_snvs_vaf{vaf}_a.vcf",
    "indel": "docs/benchmarking/snvindel/data/truth_indels_vaf{vaf}_a.vcf",
    "sv": "docs/benchmarking/sv/data/truth_svs_vaf{vaf}_a.vcf",
    "ins": "docs/benchmarking/sv/data/truth_ins_vaf{vaf}_a.vcf",
    "invdel": "docs/benchmarking/sv/data/truth_invdel_vaf{vaf}_a.vcf",
}


def get_variant_type(name):
    """Parse variant type from sample name prefix."""
    for t in ("invdel", "ins", "sv", "indel", "snv"):
        if name.startswith(t + "_"):
            return t
    return None


def parse_name(name):
    """Parse sample name into components.

    Returns (vtype, vaf_int, param_variant, replicate).
    """
    vtype = get_variant_type(name)
    if vtype is None:
        return None, None, None, None

    vaf_m = re.search(r"vaf(\d+)", name)
    if not vaf_m:
        return vtype, None, None, None
    vaf_int = int(vaf_m.group(1))

    parts_after_vaf = name.split(f"vaf{vaf_m.group(1)}_")[1]
    param_m = re.match(r"(cov\d+|fam\d+|pcr\d+|frag\d+)_(.+)", parts_after_vaf)
    if param_m:
        param_variant = param_m.group(1)
        replicate = param_m.group(2)
    else:
        param_variant = "base"
        replicate = parts_after_vaf if parts_after_vaf else "a"

    return vtype, vaf_int, param_variant, replicate


def compute_seed(vtype, vaf_int, replicate):
    """Compute deterministic seed from type, VAF, and replicate."""
    type_base = TYPE_BASES.get(vtype, 0)
    rep_offset = REP_OFFSETS.get(replicate, 0)
    return type_base + rep_offset + vaf_int


def get_params(param_variant):
    """Get parameter values for a given param sweep variant."""
    coverage = 5000.0
    family_size_mean = 4.0
    pcr_cycles = 8
    frag_mean = 167.0

    if param_variant in COVERAGES:
        coverage = COVERAGES[param_variant]
    elif param_variant in FAM_SIZES:
        family_size_mean = FAM_SIZES[param_variant]
    elif param_variant in PCR_CYCLES:
        pcr_cycles = PCR_CYCLES[param_variant]
    elif param_variant in FRAG_MEANS:
        frag_mean = FRAG_MEANS[param_variant]

    return coverage, family_size_mean, pcr_cycles, frag_mean


def generate_config(name):
    """Generate config YAML content for a sample name."""
    vtype, vaf_int, param_variant, replicate = parse_name(name)
    if vtype is None or vaf_int is None:
        return None

    vaf_fraction = vaf_int / 10000.0
    vaf_percent = vaf_fraction * 100
    purity = vaf_fraction * 2.0  # purity = 2 * VAF for diploid het

    seed = compute_seed(vtype, vaf_int, replicate)
    coverage, family_size_mean, pcr_cycles, frag_mean = get_params(param_variant)

    vaf_str = f"{vaf_int:04d}"
    mutations_vcf = MUTATION_VCFS[vtype].format(vaf=vaf_str)

    # SV types use a different reference
    if vtype in ("sv", "ins", "invdel"):
        reference = "docs/benchmarking/sv/data/ref.fa"
    else:
        reference = "docs/benchmarking/snvindel/data/ref.fa"

    upper_name = name.upper()
    comment_vaf = f"{vaf_percent:.3f}%"

    config = f"""# varforge config: {upper_name} — VAF {comment_vaf}
reference: {reference}

output:
  directory: docs/benchmarking/ml/results/sim_{name}
  fastq: true
  bam: false
  truth_vcf: true
  manifest: true

sample:
  name: {upper_name}
  read_length: 150
  coverage: {coverage:.1f}
  platform: illumina

fragment:
  model: normal
  mean: {frag_mean:.1f}
  sd: 30.0

quality:
  mean_quality: 37
  tail_decay: 0.002

tumour:
  purity: {purity:.6f}
  ploidy: 2

mutations:
  vcf: {mutations_vcf}

umi:
  length: 5
  duplex: true
  pcr_cycles: {pcr_cycles}
  family_size_mean: {family_size_mean:.1f}
  family_size_sd: 1.5

chromosomes:
  - chr1

seed: {seed}
"""
    return config


def needs_config(name):
    """Return True if this sim dir needs a config generated."""
    config_path = CONFIGS_DIR / f"{name}.yaml"
    return not config_path.exists()


def main():
    CONFIGS_DIR.mkdir(parents=True, exist_ok=True)

    # Find all sim dirs that need configs (don't have FASTQs)
    sim_dirs = sorted(d for d in RESULTS.iterdir() if d.is_dir() and d.name.startswith("sim_"))
    names_needing_config = [
        d.name[4:] for d in sim_dirs
        if not any(d.glob("*_R1.fastq.gz")) and needs_config(d.name[4:])
    ]
    # Also generate for any sim dir regardless
    all_names = [d.name[4:] for d in sim_dirs]

    generated = 0
    skipped = 0
    failed = 0

    for name in all_names:
        config_path = CONFIGS_DIR / f"{name}.yaml"
        if config_path.exists():
            skipped += 1
            continue

        config = generate_config(name)
        if config is None:
            print(f"[SKIP] Cannot generate config for {name}")
            failed += 1
            continue

        config_path.write_text(config)

        # Also write cmd.txt
        cmd_path = CONFIGS_DIR / f"{name}_cmd.txt"
        cmd_path.write_text(
            f"# Run from the repository root.\n"
            f"varforge simulate --config docs/benchmarking/ml/configs/{name}.yaml\n"
        )
        generated += 1

    print(f"Generated={generated} Skipped={skipped} Failed={failed}")
    print(f"Total configs in dir: {len(list(CONFIGS_DIR.glob('*.yaml')))}")


if __name__ == "__main__":
    main()
