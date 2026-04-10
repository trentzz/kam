# ML Models

kam includes a built-in ML scorer that re-ranks variant calls using a gradient-boosted classifier.
The scorer runs after the statistical calling stage and adds two output columns:

- `ml_prob` — model confidence (0–1) that the call is a real variant.
- `ml_filter` — PASS if `ml_prob` ≥ threshold, FAIL otherwise.

Enable it with `--ml-model`:

```bash
kam run \
  --r1 sample_R1.fastq.gz \
  --r2 sample_R2.fastq.gz \
  --targets panel_targets.fa \
  --output-dir results/ \
  --ml-model single-strand-v1
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

---

## Using a custom model

Pass a file path instead of a name to load a custom ONNX model:

```bash
kam run ... --ml-model /path/to/custom_model.onnx
```

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

Total binary overhead from all bundled models: ~656 KB.
