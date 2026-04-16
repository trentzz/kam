# ML Models

kam includes a built-in ML scorer that re-ranks variant calls using a gradient-boosted classifier.
The scorer runs after the statistical calling stage and adds two output columns:

- `ml_prob` — model confidence (0–1) that the call is a real variant.
- `ml_filter` — PASS if `ml_prob` ≥ threshold, FAIL otherwise.

Enable it with `--ml-model`. For discovery mode on Twist duplex panels, use `twist-duplex-v2`:

```bash
kam run \
  --r1 sample_R1.fastq.gz \
  --r2 sample_R2.fastq.gz \
  --targets panel_targets.fa \
  --output-dir results/ \
  --ml-model twist-duplex-v2
```

To see all available built-in models:

```bash
kam models list
```

---

## Available models

| Name | Chemistry | Training data | AUPRC | AUROC | When to use |
|---|---|---|---|---|---|
| `single-strand-v1` | Twist UMI duplex | ML3 single-strand consensus, 9,990 samples | 0.9998 | 0.9949 | Standard Twist UMI panels with single-strand consensus calling |
| `twist-duplex-v1` | Twist UMI duplex | Duplex-mode calls, 10,000 samples | 0.6062 | 0.9001 | Twist UMI panels run with duplex consensus calling enabled |
| `twist-duplex-v2` | Twist UMI duplex | Real-data retrained on 24-sample titration dataset, 49 features | 0.973 | — | Recommended for discovery mode on Twist duplex panels |

### `single-strand-v1`

**Architecture**: LightGBM (200 estimators, max depth 6, 31 leaves).

**Training set**: 9,750,546 rows from 9,990 varforge-simulated samples covering five variant classes (SNV, INDEL, INS, SV, INVDEL) at VAF 0.03–14%. Samples were drawn from a wide range of coverages (500–18,000x), fragment size distributions, and PCR cycle counts.

**Test set**: 989,225 rows from 1,000 held-out samples (same generation, no overlap with training).

**Features** (33 total, all computable from `VariantCall` at inference time):

| Feature | Description |
|---|---|
| `vaf` | Variant allele frequency |
| `nref` / `nalt` | Reference / alt molecule count |
| `ndupalt` / `nsimalt` | Duplex / simplex alt molecule count |
| `sbp` | Strand bias p-value |
| `conf` | Posterior confidence (beta-binomial) |
| `ref_len` / `alt_len` | Reference / alt allele length |
| `duplex_frac` | Fraction of alt molecules with duplex support |
| `has_duplex` | 1 if any duplex alt molecule present |
| `ci_width` | 95% CI width on VAF estimate |
| `alt_depth` | Total depth at alt position |
| `log_nalt` / `log_nref` / `log_alt_depth` / `log_vaf` | Log-transformed counts |
| `vaf_times_conf` / `vaf_times_nalt` / `nalt_over_conf` | Interaction features |
| `ci_width_rel` | CI width relative to VAF |
| `snr` | Signal-to-noise ratio (`nalt / (nref + 1)`) |
| `conf_sq` / `nalt_sq` / `vaf_sq` | Squared features |
| `ref_alt_len_ratio` | `ref_len / max(alt_len, 1)` |
| `indel_size` | `abs(alt_len - ref_len)` |
| `duplex_enrichment` | `duplex_frac / max(1 - duplex_frac, 0.01)` |
| `simplex_only_frac` | `nsimalt / max(nalt, 1)` |
| `conf_above_99` / `conf_above_999` | Binary confidence thresholds |
| `sbp_above_05` | Binary strand bias threshold |
| `variant_class_enc` | Variant class encoded as integer |

**Performance on held-out test split**:

| Metric | Value |
|---|---|
| AUPRC | 0.9998 |
| AUROC | 0.9949 |
| Precision | 0.9978 |
| Recall | 0.9966 |

### `twist-duplex-v1`

**Architecture**: LightGBM (300 estimators, max depth 6, 31 leaves).

**Training set**: 5,891,908 rows from 10,000 varforge-simulated duplex samples (535,628 positives, 9.1% positive rate after 10:1 negative subsampling). Samples span SNV, insertion, deletion, large deletion, tandem duplication, inversion, and invdel variant classes at VAF 0.03–14%, across a range of coverages, fragment sizes, and PCR cycle counts.

**Test set**: 1,541,242 rows from 1,000 held-out samples (44,852 positives, 2.9% positive rate — unsubsampled, reflecting realistic class imbalance).

**Features**: same 33 features as `single-strand-v1` (see table above), all computable from `VariantCall` at inference time.

**Performance on held-out test split**:

| Metric | Value |
|---|---|
| AUPRC | 0.6062 |
| AUROC | 0.9001 |
| Precision | 0.1263 |
| Recall | 0.7978 |

The lower AUPRC compared to `single-strand-v1` reflects the harder setting: duplex-mode calls carry greater intrinsic noise and the test set is evaluated at full class imbalance (97.1% negatives). AUROC remains high at 0.90, showing the model ranks true positives well across the score range. At the default threshold of 0.5, recall is high (0.80) at the cost of precision (0.13); lower the threshold to trade recall for precision.

### `twist-duplex-v2`

**Architecture**: LightGBM, retrained on real data with GroupKFold cross-validation (4 folds, grouped by sample ID) and 50 Optuna hyperparameter trials.

**Training set**: Real confirmed TP/FP calls from the 24-sample Twist cfDNA Pan-Cancer titration dataset (3 input masses, 8 VAF levels). Unlike v1 (trained on varforge simulations), v2 uses real sequencing data where the distinction between true somatic variants and background noise is observable.

**Features** (49 total in v2; 51 total with locus difficulty features): the 33 features from v1 plus 16 sequence-context features including substitution type, trinucleotide context, CpG status, GC content, and homopolymer run length. Two additional locus difficulty features (`dust_score` and `repeat_fraction`) are available for new model training (see below).

**Key improvement over v1**: `twist-duplex-v1` assigned ml_prob >= 0.9 to all PASS calls (both TP and FP) because it was trained on simulated data where all statistical PASS calls are real variants. `twist-duplex-v2` produces discriminating probability distributions that separate true positives from false positives. FP counts drop by 79-93% across negative controls while sensitivity is fully preserved.

**Performance**: AUPRC 0.973 on the held-out 30 ng test set. Precision >= 0.969 at 1-2% VAF. See `docs/project/experiments/03-ml-twist-duplex-v2/results/` for full evaluation.

**Default threshold**: 0.449 (optimised from training data).

---

## Using a custom model

Use `--custom-ml-model` with a file path to load an external ONNX model:

```bash
kam run ... --custom-ml-model /path/to/custom_model.onnx
```

The `--custom-ml-model` flag is mutually exclusive with `--ml-model` (which selects a built-in model by name).

The companion metadata file must exist at the same path with a `.json` extension
(e.g. `/path/to/custom_model.json`). It must have this structure:

```json
{
  "version": "3",
  "feature_names": ["vaf", "nref", ...],
  "ml_pass_threshold": 0.5,
  "variant_class_map": {"SNV": 0, "Insertion": 1, ...}
}
```

---

## Threshold

The default pass threshold is `0.5`. Calls with `ml_prob >= 0.5` receive `ml_filter = PASS`.
Lower the threshold to increase recall at the cost of precision. The threshold is read from the
model's `meta.json` and can be overridden in a future release.

---

## Model size and bundling

Built-in models are compiled into the `kam` binary using `include_bytes!`. The current binary
includes:

| Model | ONNX size |
|---|---|
| `single-strand-v1` | 656 KB |
| `twist-duplex-v1` | 668 KB |
| `twist-duplex-v2` | 2.3 MB |

Total binary overhead from all bundled models: ~3.6 MB.

---

## Locus difficulty features (features 50–51)

Two features quantify sequence complexity at the variant locus. These are available for new
model training but are not used by the current bundled models (`single-strand-v1`,
`twist-duplex-v1`, `twist-duplex-v2`). When a model does not request these features in its
metadata, they default to 0.0 in the input vector.

| Feature | Name | Description |
|---|---|---|
| 50 | `dust_score` | DUST algorithm score over a 64 bp window centred on the variant. Computed from trinucleotide frequencies. Higher values indicate more repetitive sequence. Score < 2 is complex; 2–10 is moderate; > 10 is low-complexity (homopolymers, simple repeats). |
| 51 | `repeat_fraction` | Fraction of bases in the reference window within a homopolymer run of ≥ 3 bases or a dinucleotide repeat of ≥ 6 bases. Range 0.0–1.0. A poly-A tract gives a high value; complex exonic sequence gives 0.0. |

Both features are computed from `VariantCall.ref_sequence` without external reference data.

**Why these features matter**: low-complexity regions generate excess background errors from
polymerase slippage and sequencing artefacts. A call at 0.5% VAF in a poly-T tract is far less
credible than the same call in complex exonic sequence. Without a direct complexity measure, the
model must infer this indirectly from homopolymer run length and GC content. DUST score and
repeat fraction capture the global complexity signal that those features miss.

**Backward compatibility**: existing bundled models ignore unknown feature names. The feature
extractor checks the model's metadata JSON for requested features. If `dust_score` or
`repeat_fraction` are not listed, they are not computed and default to 0.0. New models trained
with these features will list them in their `feature_names` array.
