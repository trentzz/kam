# SV Parameter Tuning

## Overview

SV detection sensitivity depends on three main parameters: k-mer size, minimum confidence
threshold, and minimum alt molecule count. The optimal values differ by SV type and VAF range.
kam provides sweep scripts to systematically explore the parameter space and identify the best
operating point for your panel and use case.

---

## Tunable parameters

### K-mer size (`-k`)

Default: 31. Range: 21 to 41 (odd values preferred).

Smaller k increases sensitivity at low VAF because shorter k-mers appear in more reads. The
tradeoff is more spurious branches in the de Bruijn graph, which increases false paths.
Larger k improves specificity but requires longer reads to span each k-mer window.

For 150 bp reads and 100 bp target windows, k values from 21 to 41 are geometrically safe
(2k does not exceed the target length).

### SV minimum confidence (`--sv-min-confidence`)

Default: 0.95. Range: 0.0 to 1.0.

The posterior probability threshold for SV PASS calls. Lower values accept more calls at the
cost of more false positives. The SV threshold is intentionally lower than the SNV threshold
(0.99) because the probability of a random error mimicking a 50+ bp structural event is
negligible.

Swept values: 0.80, 0.90, 0.95, 0.99.

### SV minimum alt molecules (`--sv-min-alt-molecules`)

Default: 1. Range: 1 to 5.

The minimum number of molecules supporting the alt allele for a PASS SV call. Setting to 1
maximises sensitivity in monitoring mode. Setting to 2 or 3 improves specificity in discovery
mode.

Swept values: 1, 2, 3.

---

## Running the k-mer sweep

The k-mer sweep re-indexes the data at each k value and re-runs the full pipeline. This is
the slower of the two sweeps because indexing is the most expensive stage.

### Script location

```
docs/benchmarking/03-sv-extended/scripts/run_kmer_sweep.sh
```

### Usage

```bash
# Set the FASTQ directory
export KAM_FASTQ_DIR=/path/to/sv_titration_fastqs

# Run the sweep
bash docs/benchmarking/03-sv-extended/scripts/run_kmer_sweep.sh
```

The script iterates over k = 21, 25, 27, 31, 35, 41 at six key VAF levels (0.1%, 0.25%,
0.5%, 1%, 2%, 5%) across all SV types and both replicates.

### Output structure

```
bigdata/benchmarking/sv_sweep/kmer_size/
├── k21/
│   ├── titration_results_discovery.tsv
│   └── titration_results_ti.tsv
├── k25/
│   └── ...
├── k31/
│   └── ...
└── k41/
    └── ...
```

Each TSV contains per-variant results with VAF, n_alt, confidence, and filter status.

---

## Running the threshold sweep

The threshold sweep re-runs only the calling stage with different confidence and molecule
thresholds. The index does not change, so this sweep is fast.

### Script location

```
docs/benchmarking/03-sv-extended/scripts/run_threshold_sweep.sh
```

### Usage

```bash
# Set the results directory from a prior benchmark run
export SV_RESULTS_DIR=/path/to/sv_baseline_results

# Run the sweep
bash docs/benchmarking/03-sv-extended/scripts/run_threshold_sweep.sh
```

The script iterates over a 4 x 3 parameter grid: sv_min_confidence in {0.80, 0.90, 0.95,
0.99} crossed with sv_min_alt_molecules in {1, 2, 3}. Each combination is evaluated on the
same six key VAF levels and both replicates.

### Output structure

```
bigdata/benchmarking/sv_sweep/threshold/
├── conf0.80_mol1/
│   ├── titration_results_discovery.tsv
│   └── titration_results_ti.tsv
├── conf0.80_mol2/
│   └── ...
├── conf0.95_mol1/
│   └── ...
└── conf0.99_mol3/
    └── ...
```

---

## Interpreting heatmap results

The scoring script produces heatmaps from sweep results:

```bash
python3 docs/benchmarking/03-sv-extended/scripts/plot_sv_new_results.py \
  --sweep-dir bigdata/benchmarking/sv_sweep/ \
  --output-dir bigdata/benchmarking/sv_sweep/figures/
```

### K-mer size vs VAF heatmap

Each cell shows the sensitivity (fraction of true positives detected) at a given k and VAF.
Cells are coloured from red (0%) through yellow (50%) to green (100%).

What to look for:

- A k value that reaches 100% sensitivity at the lowest possible VAF.
- Consistent performance across SV types. If one SV type drops at high k, that k is too large
  for its breakpoint structure.
- The default k=31 should be competitive. If k=25 is substantially better at low VAF, consider
  using it for your panel.

### Confidence vs min_alt_molecules heatmap

Each cell shows the F1 score (harmonic mean of sensitivity and specificity) at a given
threshold combination. This heatmap helps you choose operating points for discovery vs
monitoring.

What to look for:

- For monitoring mode: low confidence + min_alt_molecules=1 maximises sensitivity.
- For discovery mode: higher confidence + min_alt_molecules=2 balances sensitivity and
  specificity.
- If the F1 score is flat across a range of thresholds, the defaults are adequate.

---

## Recommended workflow

### Option 1: index once, sweep thresholds

If you are tuning for a specific panel and k=31 works well, index once at k=31 and sweep
only the confidence and molecule thresholds. This is the fast path.

```bash
# 1. Run the baseline benchmark at k=31
bash docs/benchmarking/03-sv-extended/scripts/run_sv_new_suite.sh

# 2. Sweep thresholds on the existing results
bash docs/benchmarking/03-sv-extended/scripts/run_threshold_sweep.sh

# 3. Plot heatmaps
python3 docs/benchmarking/03-sv-extended/scripts/plot_sv_new_results.py \
  --sweep-dir bigdata/benchmarking/sv_sweep/ \
  --output-dir bigdata/benchmarking/sv_sweep/figures/
```

### Option 2: full k-mer sweep

If you are exploring whether a different k improves SV detection on your chemistry, run the
full k-mer sweep first, pick the best k, then sweep thresholds at that k.

```bash
# 1. Sweep k-mer sizes
bash docs/benchmarking/03-sv-extended/scripts/run_kmer_sweep.sh

# 2. Pick the best k from the heatmap
# 3. Re-run threshold sweep at the chosen k
KAM_KMER_SIZE=25 bash docs/benchmarking/03-sv-extended/scripts/run_threshold_sweep.sh
```

---

## Notes

- The k-mer sweep takes 6 to 12 hours depending on hardware because it re-indexes at each k.
- The threshold sweep takes minutes because it reuses existing index and path data.
- All sweep scripts write to `bigdata/` and are not committed to the repository.
- Use the same FASTQ data for all sweeps. Mixing data sources invalidates comparisons.
