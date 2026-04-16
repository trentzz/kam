# ML Locus Difficulty Features

## Status

In Progress.

## Summary

Two new ML features (features 50 and 51) quantifying sequence complexity at
the variant locus. DUST score detects low-complexity windows via trinucleotide
frequency analysis. Repeat fraction measures the fraction of bases in
homopolymer or dinucleotide repeats.

## Problem

The v2 ML model had no direct measure of global sequence complexity beyond GC
content and homopolymer run length. Low-complexity regions generate excess
background errors that are difficult to distinguish from true variants without
locus-specific difficulty signals. A call at 0.5% VAF in a poly-T tract is
far less credible than the same call in a complex exonic region, but the model
had no feature to capture this difference.

## Solution

Two new features, both computable from the `VariantCall.ref_sequence` field
without any external reference data:

### Feature 50: `dust_score`

DUST algorithm applied to a sliding 64 bp window centred on the variant
position. The score is derived from trinucleotide frequency analysis within
the window: for each trinucleotide, count occurrences, compute
`count * (count - 1) / 2`, and sum over all trinucleotides. Divide by the
window length minus 2. Higher scores indicate more repetitive sequence.

Interpretation:
- Score < 2: complex, unique sequence
- Score 2–10: moderate complexity
- Score > 10: low-complexity region (homopolymer runs, simple repeats)

### Feature 51: `repeat_fraction`

Fraction of bases in the reference window that fall within a homopolymer run
of 3 or more bases, or a dinucleotide repeat of 6 or more bases (3 or more
repeat units). Value ranges from 0.0 to 1.0.

Examples:
- `AAAAA` in a 20 bp window: 5/20 = 0.25 contribution from that run
- `ATATAT` in a 20 bp window: 6/20 = 0.30 contribution from that run
- Complex sequence with no repeats: 0.0

### Backward compatibility

Existing bundled models (`single-strand-v1`, `twist-duplex-v1`,
`twist-duplex-v2`) specify their feature names in metadata JSON. Unknown
feature names are ignored by the feature extractor and default to 0.0 in the
ONNX input vector. New models trained with these features will list them in
their metadata; old models will not request them.

## Tests

- Poly-A input gives DUST score > 10
- Complex exonic sequence gives DUST score < 2
- Poly-A input gives repeat fraction > 0.5
- Short sequence with no repeats gives repeat fraction 0.0
- Dinucleotide repeat (ATATAT...) gives repeat fraction > 0.5
- Feature extraction with an old model metadata file produces 0.0 for both
  new features (backward compatibility)

## Relevant Files

- `kam-call/src/ml/features.rs` — feature extraction logic
- `kam-call/src/ml/model.rs` — ONNX feature vector construction
- `kam-call/src/ml/metadata.rs` — model metadata parsing
