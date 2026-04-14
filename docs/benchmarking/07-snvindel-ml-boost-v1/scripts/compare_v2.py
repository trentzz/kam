#!/usr/bin/env python3
"""Compare discovery, disc+ML_v1, disc+ML_v2, and TI mode results.

Reads:
  results/titration_2mreads_disc.tsv
  results/titration_2mreads_disc_ml.tsv
  results/titration_2mreads_disc_ml_v2.tsv
  results/titration_2mreads_ti_corrected.tsv

Prints a table and writes a CSV summary for each metric across conditions.

Run from the repository root:
    python3 docs/benchmarking/07-snvindel-ml-boost-v1/scripts/compare_v2.py
"""

from pathlib import Path
import csv

RESULTS = Path("docs/benchmarking/07-snvindel-ml-boost-v1/results")

DISC_TSV    = RESULTS / "titration_2mreads_disc.tsv"
V1_TSV      = RESULTS / "titration_2mreads_disc_ml.tsv"
V2_TSV      = RESULTS / "titration_2mreads_disc_ml_v2.tsv"
TI_TSV      = RESULTS / "titration_2mreads_ti_corrected.tsv"


def load_tsv(path):
    rows = {}
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            ng  = row["ng"]
            vaf = row["vaf"]
            rows[(ng, vaf)] = row
    return rows


def fmt(v, digits=3):
    try:
        return f"{float(v):.{digits}f}"
    except (ValueError, TypeError):
        return str(v)


def main():
    disc = load_tsv(DISC_TSV)
    v1   = load_tsv(V1_TSV)
    v2   = load_tsv(V2_TSV)
    ti   = load_tsv(TI_TSV)

    vafs = ["0.0", "0.001", "0.01", "0.1", "0.25", "0.5", "1.0", "2.0"]
    ngs  = ["5ng", "15ng", "30ng"]

    print("=" * 120)
    print("PRECISION comparison  (goal: disc+ML_v2 approaches TI precision=1.000 with minimal sensitivity loss)")
    print("=" * 120)
    print(f"{'ng':>6} {'vaf%':>6} | {'Disc prec':>10} {'v1 prec':>10} {'v2 prec':>10} {'TI prec':>10} "
          f"| {'Disc FP':>8} {'v1 FP':>6} {'v2 FP':>6} {'TI FP':>6} "
          f"| {'Disc sens':>10} {'v2 ml_sens':>12}")
    print("-" * 120)

    rows_out = []

    for ng in ngs:
        for vaf in vafs:
            key = (ng, vaf)
            d  = disc.get(key, {})
            r1 = v1.get(key, {})
            r2 = v2.get(key, {})
            t  = ti.get(key, {})

            d_prec  = d.get("precision", "NA")
            v1_prec = r1.get("precision", "NA")
            # For v2: use ml_precision (precision after ML filter)
            v2_prec = r2.get("ml_precision", "NA")
            ti_prec = t.get("precision", "NA")

            d_fp    = d.get("fp", "NA")
            v1_fp   = r1.get("fp", "NA")
            v2_fp   = r2.get("ml_fp", "NA")
            ti_fp   = t.get("fp", "NA")

            d_sens   = d.get("sensitivity", "NA")
            v2_sens  = r2.get("ml_sensitivity", "NA")

            vaf_pct = f"{float(vaf)*100:.3g}" if vaf != "NA" else vaf

            print(f"{ng:>6} {vaf_pct:>6} | {fmt(d_prec):>10} {fmt(v1_prec):>10} {fmt(v2_prec):>10} {fmt(ti_prec):>10} "
                  f"| {fmt(d_fp, 0):>8} {fmt(v1_fp, 0):>6} {fmt(v2_fp, 0):>6} {fmt(ti_fp, 0):>6} "
                  f"| {fmt(d_sens):>10} {fmt(v2_sens):>12}")

            rows_out.append({
                "ng": ng, "vaf": vaf,
                "disc_precision": d_prec, "v1_precision": v1_prec, "v2_precision": v2_prec, "ti_precision": ti_prec,
                "disc_fp": d_fp, "v1_fp": v1_fp, "v2_fp": v2_fp, "ti_fp": ti_fp,
                "disc_sensitivity": d_sens, "v2_ml_sensitivity": v2_sens,
                "ti_sensitivity": t.get("sensitivity", "NA"),
            })

    print()
    print("=" * 120)
    print("KEY NEGATIVE CONTROLS (vaf=0): FP counts")
    print("=" * 120)
    for ng in ngs:
        key = (ng, "0.0")
        d = disc.get(key, {})
        r2 = v2.get(key, {})
        t = ti.get(key, {})
        print(f"  {ng}: disc FP={d.get('fp','NA')}  v2 ML FP={r2.get('ml_fp','NA')}  TI FP={t.get('fp','NA')}")

    print()
    print("=" * 120)
    print("SENSITIVITY RETENTION at key VAF levels (v2 ml_sens / disc_sens)")
    print("=" * 120)
    for ng in ngs:
        for vaf in ["1.0", "2.0"]:
            key = (ng, vaf)
            d = disc.get(key, {})
            r2 = v2.get(key, {})
            try:
                ratio = float(r2.get("ml_sensitivity", 0)) / float(d.get("sensitivity", 1))
                print(f"  {ng} {vaf}% VAF: disc_sens={fmt(d.get('sensitivity','NA'))}  v2_ml_sens={fmt(r2.get('ml_sensitivity','NA'))}  retention={ratio:.1%}")
            except (ValueError, ZeroDivisionError):
                print(f"  {ng} {vaf}% VAF: NA")

    # Write CSV
    out_path = RESULTS / "comparison_disc_v1_v2_ti.csv"
    if rows_out:
        with open(out_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=list(rows_out[0].keys()))
            writer.writeheader()
            writer.writerows(rows_out)
        print(f"\nWrote: {out_path}")


if __name__ == "__main__":
    main()
