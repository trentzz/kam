# SV Parameter Tuning Guide

**Status**: done
**Epic**: SV-SWEEP

## Summary

Empirical parameter sweep across k-mer size, confidence thresholds, minimum
alt molecule counts, and SV size to find optimal operating points for each
SV type.

## Parameters Swept

- k-mer size: 21, 25, 27, 31, 35, 41
- sv_min_confidence: 0.80, 0.90, 0.95, 0.99
- sv_min_alt_molecules: 1, 2, 3
- SV size: 20, 50, 100, 200, 500 bp

## Deliverables

- Heatmap: k-mer size vs VAF, coloured by sensitivity, per SV type
- Heatmap: confidence threshold vs min_alt_molecules
- Line plot: SV size vs sensitivity per type
- Recommended configuration per SV type
