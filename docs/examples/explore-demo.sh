#!/usr/bin/env bash
# Explore demo: generate synthetic Twist UMI duplex data with varforge,
# assemble molecules with kam, then open the REPL to inspect results.
#
# Usage: bash docs/examples/explore-demo.sh
# Requires: varforge >= 0.2.0, kam (built in release mode)

set -euo pipefail

DEMO_DIR=~/tmp/kam-explore-demo
CONFIG=docs/benchmarking/01-snvindel/configs/snv_vaf010.yaml

echo "==> Generating synthetic data with varforge..."
mkdir -p "$DEMO_DIR"
varforge simulate \
    --config "$CONFIG" \
    --output-dir "$DEMO_DIR/sim" \
    --preset umi \
    --coverage 2000 \
    --seed 42

R1=$(ls "$DEMO_DIR"/sim/*_R1*.fastq.gz | head -1)
R2=$(ls "$DEMO_DIR"/sim/*_R2*.fastq.gz | head -1)

echo "==> Assembling molecules..."
kam assemble \
    --r1 "$R1" \
    --r2 "$R2" \
    --output "$DEMO_DIR/molecules.bin" \
    --chemistry twist-umi-duplex

echo ""
echo "==> Opening REPL. Try: summary, stats, head 20, filter n_duplex_reads > 4, quit"
echo ""
kam explore "$DEMO_DIR/molecules.bin"
