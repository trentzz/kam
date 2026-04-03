#!/usr/bin/env bash
# Run the full kam SV detection pipeline for all VAF levels.
#
# Runs two passes per VAF level:
#   1. Discovery mode  — no --target-variants; all variants reported.
#   2. Monitoring mode — --target-variants derived VCF; FPs relabelled NotTargeted.
#
# The monitoring truth VCF is generated from the discovery PASS SV calls using
# make_monitoring_vcf.py, which applies the same variant-key extraction logic
# as kam's internal apply_target_filter. This guarantees exact key matching
# regardless of VCF representation differences between the simulation truth VCF
# and kam's internal normalisation.
#
# Usage: bash run_kam_sv.sh
#
# Requires: kam binary in target/release/kam, Python 3
#
# Outputs:
#   docs/benchmarking/02-sv-core/results/kam_vaf{tag}/            — discovery
#   docs/benchmarking/02-sv-core/results/kam_vaf{tag}_monitoring/ — monitoring

set -euo pipefail

KAM=./target/release/kam
EXPDIR=docs/benchmarking/sv
DATADIR=$EXPDIR/data
RESDIR=$EXPDIR/results
SCRIPTS=$EXPDIR/scripts

TARGETS=$DATADIR/sv_targets.fa
JUNCTIONS=$DATADIR/sv_junctions.fa

for TAG in 005 010 020 050; do
    SIMDIR=$RESDIR/sim_vaf${TAG}
    OUTDIR=$RESDIR/kam_vaf${TAG}
    MONDIR=$RESDIR/kam_vaf${TAG}_monitoring
    MONTRUTH=$OUTDIR/monitoring_truth.vcf

    mkdir -p "$OUTDIR" "$MONDIR"

    R1=$(ls "$SIMDIR"/*_R1.fastq.gz)
    R2=$(ls "$SIMDIR"/*_R2.fastq.gz)

    echo "=== VAF ${TAG} — discovery ==="
    $KAM run \
        --r1 "$R1" \
        --r2 "$R2" \
        --targets "$TARGETS" \
        --sv-junctions "$JUNCTIONS" \
        --kmer-size 31 \
        --output-format tsv \
        --output-dir "$OUTDIR" \
        --strand-bias-threshold 0 \
        2>&1 | grep -E "^\[|ERROR|molecules|on.target" || true
    echo "  Done: $OUTDIR/"

    # Generate monitoring truth VCF from discovery PASS SV calls.
    python3 "$SCRIPTS/make_monitoring_vcf.py" "$OUTDIR/variants.tsv" "$MONTRUTH"

    echo "=== VAF ${TAG} — monitoring ==="
    $KAM run \
        --r1 "$R1" \
        --r2 "$R2" \
        --targets "$TARGETS" \
        --sv-junctions "$JUNCTIONS" \
        --kmer-size 31 \
        --output-format tsv \
        --output-dir "$MONDIR" \
        --strand-bias-threshold 0 \
        --target-variants "$MONTRUTH" \
        2>&1 | grep -E "^\[|ERROR|molecules|on.target" || true
    echo "  Done: $MONDIR/"
done

echo ""
echo "All VAFs complete."
