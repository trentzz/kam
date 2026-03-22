#!/usr/bin/env bash
# Run the full kam SV detection pipeline for all VAF levels.
#
# Usage: bash run_kam_sv.sh
#
# Requires: kam binary in target/release/kam
#           pyfaidx installed (for make_sv_targets.py)
#
# Outputs: docs/benchmarking/sv/results/kam_vaf{tag}/

set -euo pipefail

KAM=./target/release/kam
EXPDIR=docs/benchmarking/sv
DATADIR=$EXPDIR/data
RESDIR=$EXPDIR/results

TARGETS=$DATADIR/sv_targets.fa
JUNCTIONS=$DATADIR/sv_junctions.fa

for TAG in 005 010 020 050; do
    SIMDIR=$RESDIR/sim_vaf${TAG}
    OUTDIR=$RESDIR/kam_vaf${TAG}
    mkdir -p "$OUTDIR"

    R1=$(ls "$SIMDIR"/*_R1.fastq.gz)
    R2=$(ls "$SIMDIR"/*_R2.fastq.gz)

    echo "=== VAF ${TAG} ==="
    echo "  Running..."
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
done

echo ""
echo "All VAFs complete."
