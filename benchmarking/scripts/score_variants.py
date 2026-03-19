#!/usr/bin/env python3
"""
Score kam VCF output against a truth VCF.

Matching is by (REF, ALT) only. kam is alignment-free and reports POS=1 for
all variants, so position-based matching is not applicable.

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
from typing import Set, Tuple


VariantKey = Tuple[str, str]  # (REF, ALT)


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


def load_vcf_variants(path: str) -> Set[VariantKey]:
    """
    Parse a VCF file and return a set of (REF, ALT) tuples.

    Skips header lines (starting with #). Handles multi-allelic ALT fields by
    splitting on comma and creating one tuple per allele.
    """
    variants: Set[VariantKey] = set()
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
            for alt in alt_field.split(","):
                alt = alt.strip()
                if alt and alt != ".":
                    variants.add((ref, alt))
    return variants


def compute_metrics(
    truth: Set[VariantKey],
    called: Set[VariantKey],
) -> dict:
    """
    Compute TP, FP, FN and derived metrics.

    TN is defined as the number of truth variants not called, mirroring the
    convention used when evaluating against a known truth set of fixed size.
    """
    tp = len(truth & called)
    fp = len(called - truth)
    fn = len(truth - called)
    # TN: truth variants that are genuinely absent and not called.
    # For a closed truth set, TN = |truth| (all truth variants are either
    # called (TP) or missed (FN); TN counts those correctly left uncalled).
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
    "tp",
    "fp",
    "fn",
    "tn",
    "sensitivity",
    "specificity",
    "precision",
    "f1",
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

    truth = load_vcf_variants(args.truth)
    called = load_vcf_variants(args.called)

    metrics = compute_metrics(truth, called)

    row = {
        "sample_name": args.sample_name,
        "vaf": args.vaf,
        "dna_input": args.dna_input,
        **metrics,
    }

    append_results(args.output, row)

    print(
        f"[score_variants] {args.sample_name}: "
        f"TP={metrics['tp']} FP={metrics['fp']} FN={metrics['fn']} "
        f"sensitivity={metrics['sensitivity']:.4f} "
        f"precision={metrics['precision']:.4f} "
        f"F1={metrics['f1']:.4f}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
