#!/usr/bin/env python3
"""
Convert the titration truth TSV to a standard VCF 4.2 file.

Input TSV columns: Number, Chromosome, Position, Ref, Alt, Type
Output: benchmarking/scripts/truth_variants.vcf
"""

import argparse
import csv
import os
import sys
from datetime import date


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Convert titration truth TSV to VCF 4.2"
    )
    parser.add_argument(
        "--input",
        default=os.path.join(
            os.path.dirname(__file__),
            "../../..",
            os.environ.get("KAM_PROBES_TSV", "titration.probes.QC.pass.tsv"),
        ),
        help="Path to input TSV file",
    )
    parser.add_argument(
        "--output",
        default=os.path.join(os.path.dirname(__file__), "truth_variants.vcf"),
        help="Path to output VCF file",
    )
    return parser.parse_args()


def convert(input_path: str, output_path: str) -> int:
    """Read the TSV and write a valid VCF 4.2 file. Returns the variant count."""
    variants = []

    with open(input_path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        # Strip surrounding quotes from header names produced by R-style quoting.
        reader.fieldnames = [f.strip('"') for f in reader.fieldnames]

        for row in reader:
            # Values may also be quoted (e.g. "G").
            chrom = row["Chromosome"].strip('"')
            pos = int(row["Position"].strip('"'))
            ref = row["Ref"].strip('"')
            alt = row["Alt"].strip('"')
            variant_type = row["Type"].strip('"').upper()

            # Normalise type to SNP or INDEL for the INFO field.
            if variant_type not in ("SNP", "INDEL"):
                variant_type = "INDEL" if len(ref) != len(alt) else "SNP"

            variants.append((chrom, pos, ref, alt, variant_type))

    # Sort by chromosome then position for a well-formed VCF.
    def sort_key(v):
        chrom = v[0]
        # Sort chrN numerically where possible.
        num = chrom.replace("chr", "")
        try:
            return (0, int(num), v[1])
        except ValueError:
            # chrX, chrY, chrM, etc.
            return (1, num, v[1])

    variants.sort(key=sort_key)

    today = date.today().strftime("%Y%m%d")

    with open(output_path, "w") as out:
        out.write("##fileformat=VCFv4.2\n")
        out.write(f"##fileDate={today}\n")
        out.write("##reference=GRCh38\n")
        out.write('##INFO=<ID=TYPE,Number=1,Type=String,Description="Variant type: SNP or INDEL">\n')
        out.write(
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        )
        for chrom, pos, ref, alt, variant_type in variants:
            out.write(
                f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\tTYPE={variant_type}\n"
            )

    return len(variants)


def main() -> None:
    args = parse_args()
    n = convert(args.input, args.output)
    print(f"Wrote {n} variants to {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
