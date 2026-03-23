#!/usr/bin/env python3
"""
Score kam VCF output against a truth VCF.

Matching is by (REF, ALT) only. kam is alignment-free and reports POS=1 for
all variants, so position-based matching is not applicable.

Variant classification:
- SNV: REF and ALT are both single bases (len == 1).
- Indel: anything else (insertion, deletion, MNV).

The truth VCF uses INFO/TYPE=SNP or TYPE=INDEL when present. Called VCFs may
not have this tag; length-based classification is used as a fallback.

Usage:
    python score_variants.py \\
        --truth truth_variants.vcf \\
        --called called.vcf \\
        --output results.tsv \\
        --sample-name Sample_15ng_VAF_1pc \\
        --vaf 1.0 \\
        --dna-input 15
"""

import argparse
import csv
import os
import sys
from typing import Dict, Set, Tuple


VariantKey = Tuple[str, str]  # (REF, ALT)

# Variant type labels.
TYPE_SNV = "snv"
TYPE_INDEL = "indel"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Score kam variant calls against a truth VCF"
    )
    parser.add_argument("--truth", required=True, help="Path to truth VCF")
    parser.add_argument("--called", required=True, help="Path to kam output VCF")
    parser.add_argument("--output", required=True, help="Path to results TSV")
    parser.add_argument(
        "--sample-name",
        default="unknown",
        help="Sample name to write in the results row",
    )
    parser.add_argument(
        "--vaf",
        type=float,
        default=float("nan"),
        help="Expected VAF for this sample",
    )
    parser.add_argument(
        "--dna-input",
        type=float,
        default=float("nan"),
        help="DNA input in ng for this sample",
    )
    return parser.parse_args()


def classify_variant(ref: str, alt: str, info: str) -> str:
    """
    Classify a variant as TYPE_SNV or TYPE_INDEL.

    Uses the INFO/TYPE tag when present (TYPE=SNP → snv, TYPE=INDEL → indel).
    Falls back to length-based classification: single-base REF and single-base
    ALT is an SNV; everything else is an indel.
    """
    for field in info.split(";"):
        if field.startswith("TYPE="):
            val = field[5:].upper()
            if val == "SNP":
                return TYPE_SNV
            if val == "INDEL":
                return TYPE_INDEL
    # Fallback: length-based.
    if len(ref) == 1 and len(alt) == 1:
        return TYPE_SNV
    return TYPE_INDEL


def load_vcf_variants(
    path: str,
) -> Tuple[Set[VariantKey], Dict[str, Set[VariantKey]]]:
    """
    Parse a VCF file and return all variants plus a per-type breakdown.

    Returns:
        all_variants: set of (REF, ALT) tuples
        by_type: dict mapping TYPE_SNV / TYPE_INDEL to sets of (REF, ALT)

    Skips header lines (starting with #). Handles multi-allelic ALT fields by
    splitting on comma and creating one tuple per allele.
    """
    all_variants: Set[VariantKey] = set()
    by_type: Dict[str, Set[VariantKey]] = {TYPE_SNV: set(), TYPE_INDEL: set()}

    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("#") or not line:
                continue
            fields = line.split("\t")
            if len(fields) < 5:
                continue
            ref = fields[3]
            alt_field = fields[4]
            info = fields[7] if len(fields) > 7 else ""
            for alt in alt_field.split(","):
                alt = alt.strip()
                if not alt or alt == ".":
                    continue
                key: VariantKey = (ref, alt)
                all_variants.add(key)
                vtype = classify_variant(ref, alt, info)
                by_type[vtype].add(key)

    return all_variants, by_type


def compute_metrics(
    truth: Set[VariantKey],
    called: Set[VariantKey],
) -> dict:
    """
    Compute TP, FP, FN and derived metrics for a given set of variants.
    """
    tp = len(truth & called)
    fp = len(called - truth)
    fn = len(truth - called)
    tn = len(truth)

    sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    specificity = tn / (tn + fp) if (tn + fp) > 0 else 0.0
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    f1 = (
        2 * precision * sensitivity / (precision + sensitivity)
        if (precision + sensitivity) > 0
        else 0.0
    )

    return {
        "tp": tp,
        "fp": fp,
        "fn": fn,
        "tn": tn,
        "sensitivity": sensitivity,
        "specificity": specificity,
        "precision": precision,
        "f1": f1,
    }


COLUMNS = [
    "sample_name",
    "vaf",
    "dna_input",
    # Overall metrics.
    "tp",
    "fp",
    "fn",
    "tn",
    "sensitivity",
    "specificity",
    "precision",
    "f1",
    # SNV-only metrics.
    "snv_tp",
    "snv_fp",
    "snv_fn",
    "snv_sensitivity",
    "snv_precision",
    "snv_f1",
    # Indel-only metrics.
    "indel_tp",
    "indel_fp",
    "indel_fn",
    "indel_sensitivity",
    "indel_precision",
    "indel_f1",
]


def append_results(output_path: str, row: dict) -> None:
    """
    Append one result row to a TSV. Write the header if the file does not yet
    exist or is empty.
    """
    write_header = not os.path.exists(output_path) or os.path.getsize(output_path) == 0

    with open(output_path, "a", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=COLUMNS, delimiter="\t")
        if write_header:
            writer.writeheader()
        writer.writerow(row)


def main() -> None:
    args = parse_args()

    truth_all, truth_by_type = load_vcf_variants(args.truth)
    called_all, called_by_type = load_vcf_variants(args.called)

    overall = compute_metrics(truth_all, called_all)
    snv = compute_metrics(truth_by_type[TYPE_SNV], called_by_type[TYPE_SNV])
    indel = compute_metrics(truth_by_type[TYPE_INDEL], called_by_type[TYPE_INDEL])

    row = {
        "sample_name": args.sample_name,
        "vaf": args.vaf,
        "dna_input": args.dna_input,
        # Overall.
        **overall,
        # SNV.
        "snv_tp": snv["tp"],
        "snv_fp": snv["fp"],
        "snv_fn": snv["fn"],
        "snv_sensitivity": snv["sensitivity"],
        "snv_precision": snv["precision"],
        "snv_f1": snv["f1"],
        # Indel.
        "indel_tp": indel["tp"],
        "indel_fp": indel["fp"],
        "indel_fn": indel["fn"],
        "indel_sensitivity": indel["sensitivity"],
        "indel_precision": indel["precision"],
        "indel_f1": indel["f1"],
    }

    append_results(args.output, row)

    print(
        f"[score_variants] {args.sample_name}: "
        f"TP={overall['tp']} FP={overall['fp']} FN={overall['fn']} "
        f"sensitivity={overall['sensitivity']:.4f} "
        f"precision={overall['precision']:.4f} "
        f"F1={overall['f1']:.4f} | "
        f"SNV sens={snv['sensitivity']:.4f} "
        f"Indel sens={indel['sensitivity']:.4f}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
