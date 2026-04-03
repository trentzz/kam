# ML-001: Build training dataset from benchmark calls

**Epic**: ML-BOOST (overallplans/ML-BOOST.md)
**Priority**: high
**Depends on**: none
**Status**: todo

## Goal

Write a Python script that walks the SNV/indel and SV benchmark sample
directories, joins each call TSV against the truth VCF, labels calls as true
positives or false positives, derives engineered features, and emits a single
`training_data.parquet` (or CSV) ready for model training.

## Steps

1. Walk `docs/benchmarking/snvindel/samples/` and
   `docs/benchmarking/sv/samples/`.
2. For each sample, load `discovery.tsv` and `truth.tsv`.
3. Match calls to truth by chrom+pos+ref+alt (exact match). Label matched
   calls `label=1`, all others `label=0`.
4. Derive features:
   - `duplex_frac = ndupalt / nalt` (0 if nalt=0)
   - `has_duplex = ndupalt > 0`
   - `ci_width = vaf_hi - vaf_lo`
   - `ref_len = len(ref)`
   - `alt_len = len(alt)`
   - `variant_class`: SNV if ref_len==1 and alt_len==1, indel if one is >1
     and max<50, SV otherwise
   - `sample_id` (for group k-fold)
   - `dataset_type` (snvindel or sv)
5. Write to `docs/benchmarking/ml/training_data.csv`.
6. Print basic label stats: total calls, true positives, false positives,
   class balance.

## Notes

- Truth matching must be exact on chrom+pos+ref+alt. Do not fuzzy-match.
- Keep both `discovery.tsv` and `tumour_informed.tsv` — tag which mode each
  row came from with a `mode` column.
- Write the script to `scripts/ml/build_training_data.py`.

## Success criteria

- [ ] Script exists at `scripts/ml/build_training_data.py`
- [ ] `docs/benchmarking/ml/training_data.csv` exists with correct schema
- [ ] Label stats printed: sanity-check that TP rate is plausible (<<10%)
- [ ] No crash on any sample directory
