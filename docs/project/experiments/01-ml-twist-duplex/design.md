# ML Dataset for Twist Duplex Chemistry

## Motivation

The existing ML3 dataset (`docs/benchmarking/ml/`) treats Twist chemistry as one of several possible configurations. It mixes samples with variable UMI lengths, chemistry types, and inline flags. A classifier trained on that data must learn both "what does a real variant look like" and "what does Twist chemistry look like at the same time."

The Twist duplex dataset fixes all chemistry parameters and varies only the biological and sequencing noise parameters. This focuses the model entirely on the variant-signal problem. It also enables features that are only meaningful for duplex UMI data: duplex-vs-simplex VAF agreement, strand-level balance within the simplex pool, and variant-specific duplex counts. These features are meaningless or unavailable for non-duplex chemistries, so they were not included in the generic model.

The practical case for a chemistry-specific model is that deployed samples will always be Twist duplex. A model that was also trained on non-duplex data carries irrelevant variance and may underfit the duplex-specific signal structure.

## Dataset Design

### Variant types and counts

| Type | Train | Test | Notes |
|---|---|---|---|
| snv | 3,500 | 350 | Single nucleotide substitution |
| indel_short | 2,000 | 200 | 1-5bp insertion or deletion |
| indel_medium | 1,000 | 100 | 6-20bp insertion or deletion |
| sv_del | 500 | 50 | Small deletion (20-50bp) |
| sv_dup | 300 | 30 | Tandem duplication (50-100bp) |
| sv_inv | 300 | 30 | Inversion (50-100bp) |
| sv_large_del | 400 | 40 | Large deletion (100-200bp) |
| ins | 1,000 | 100 | Large insertion (70-100bp) |
| invdel | 500 | 50 | Inversion-deletion (80-130bp) |
| **Total** | **10,000** | **1,000** | |

SNV and short indel samples are overrepresented because they are the most clinically relevant variant types and the hardest to distinguish from sequencing error at low VAF.

### VAF range

VAF is sampled log-uniformly between 0.01% and 5% (i.e. `np.exp(np.random.uniform(np.log(0.0001), np.log(0.05)))`). This distribution places more samples at very low VAF where the model needs the most signal to generalise. The purity sent to varforge is `min(VAF * 2, 1.0)` for a diploid tumour.

### Parameter distributions

All samples use fixed Twist chemistry:

```yaml
umi:
  length: 5
  duplex: true
  inline: true
  family_size_sd: 1.0  # fixed
```

Varied parameters (per sample, seeded deterministically from index):

| Parameter | Distribution | Range |
|---|---|---|
| coverage | log-uniform | 1,000–15,000 |
| family_size_mean | uniform | 2.0–8.0 |
| pcr_cycles | uniform int | 5–14 |
| fragment_mean | uniform | 130–280 bp |
| fragment_sd | uniform | 20–50 bp |
| mean_quality | uniform int | 30–40 |

Each sample has 1-3 variants placed at non-overlapping positions with a minimum spacing of 150 bp (300 bp for SVs).

The deterministic seeding formula is `seed = sample_index * 7 + 42`, which ensures consistent regeneration while keeping seeds well-separated across samples.

## New Features from Updated kam-call Output

The updated kam TSV adds 7 columns beyond the original 14:

| Column | Description |
|---|---|
| `n_simplex_fwd_alt` | Simplex alt reads on the forward strand |
| `n_simplex_rev_alt` | Simplex alt reads on the reverse strand |
| `n_duplex_ref` | Duplex molecules supporting reference |
| `n_simplex_ref` | Simplex molecules supporting reference |
| `mean_alt_error_prob` | Mean base error probability at the alt position |
| `min_variant_specific_duplex` | Minimum variant-specific duplex count across replicates |
| `mean_variant_specific_molecules` | Mean variant-specific molecule count |

These columns enable a set of duplex-specific derived features that have no equivalent in the generic ML3 feature set.

## Feature Engineering

### Duplex-specific features

The most informative new features exploit the split between duplex and simplex evidence:

- `duplex_vaf` and `simplex_vaf`: VAF estimated from duplex-only and simplex-only counts. A real variant should show concordant duplex and simplex VAFs. A sequencing error is more likely to appear predominantly in simplex reads.
- `duplex_simplex_vaf_delta`: The signed difference. Systematic positive values (duplex VAF exceeds simplex) suggest a genuine variant where duplex consensus suppresses background error. Systematic negative values suggest a strand-specific artefact.
- `vaf_duplex_agreement`: How far the reported VAF deviates from the duplex-only estimate. Calibration quality signal.
- `strand_balance` and `strand_asymmetry`: Derived from `n_simplex_fwd_alt` and `n_simplex_rev_alt`. Strand bias is a classic artefact signature, but this dataset also includes the background expected balance under null conditions, so the model can learn the full distribution rather than just a threshold.
- `variant_specific_duplex_frac`: Fraction of duplex alt support that is variant-specific (not attributable to background). High values are a strong positive signal.
- `qual_confidence_interaction` and `qual_vaf_interaction`: Products of error probability with confidence and VAF. These capture the joint regime where both quality and allele frequency support a call.

### Log transforms

VAF and all count features are log-transformed because their distributions span several orders of magnitude and gradient boosters handle monotone transformations well when the raw scale is very skewed.

### Binned features

- `depth_bucket`: Quintile bin of total depth. Captures non-linear depth effects without assuming a functional form.
- `sbp_cat`: Strand bias p-value discretised into low/medium/high. Provides a coarse but robust signal for the model to anchor on.
- `conf_above_99`, `conf_above_999`: Binary flags for high-confidence calls. Allow the model to learn threshold effects that are already semantically meaningful in kam's scoring.

### Parameter features

Coverage, family size, PCR cycles, and fragment mean are included from `params.json`. At inference time, these values are not known from the FASTQ alone. They are therefore treated as nuisance covariates during training (the model should learn to be robust to them) rather than primary predictors.

## Model Design

### LightGBM and XGBoost

Both models are trained because they have complementary strengths on tabular data. LightGBM tends to be faster and handles imbalanced data better with `is_unbalance=True`. XGBoost with `scale_pos_weight` handles class imbalance more explicitly. The ensemble averages their calibrated probabilities with equal weight. Equal weighting is a principled starting point; a stacked meta-learner could be added once the dataset is large enough to hold out a validation set for it.

### Optuna HPO

200 Optuna TPE trials per model optimise AUPRC using 3-fold GroupKFold CV (grouped by sample ID) for speed during the search. The full 5-fold CV runs only once with the best parameters, to measure generalisation variance without inflating the HPO budget.

AUPRC is the optimisation target rather than AUROC because the dataset is class-imbalanced (most positions are not true variants) and AUPRC is more sensitive to performance at the positive end of the score distribution.

### Calibration

Both models are calibrated with `CalibratedClassifierCV(method='isotonic', cv='prefit')` on the full training set. Isotonic regression is preferred over Platt scaling because the score distributions from gradient boosters are not sigmoid-shaped. The calibrated probabilities can be interpreted as posterior variant probabilities for downstream decision-making.

### ONNX export

The XGBoost model is exported to ONNX for deployment. The LightGBM model is saved in its native text format (which can be loaded in any environment with LightGBM installed). The ensemble is saved as a Python pickle for evaluation and comparison purposes only; the ONNX model is intended for production use.

## Placeholder Results

*Note: This section will be filled after training runs complete.*

| Model | AUROC | AUPRC | Sens @ 1% FPR | Best F1 |
|---|---|---|---|---|
| LightGBM | — | — | — | — |
| XGBoost | — | — | — | — |
| Ensemble | — | — | — | — |

### Per-VAF-bin AUPRC (Ensemble)

| VAF Bin | AUPRC | N samples |
|---|---|---|
| Very low (<0.1%) | — | — |
| Low (0.1–0.5%) | — | — |
| Medium (0.5–2%) | — | — |
| High (>2%) | — | — |

## References

- Dataset generation: `scripts/ml/generate_twist_duplex_configs.py`
- Pipeline runner: `scripts/ml/run_twist_duplex_pipeline.py`
- Feature extraction: `scripts/ml/build_twist_duplex_features.py`
- Model training: `scripts/ml/train_twist_duplex.py`
- Nextcloud upload: `scripts/ml/upload_twist_duplex.sh`
- Experiment README: `docs/project/experiments/01-ml-twist-duplex/README.md`
- kam-call output columns: `kam-call/src/output.rs`
- Twist chemistry details: `docs/research/twist_umi_chemistry.md`
- Statistical calling models: `docs/research/statistical_calling_models.md`
