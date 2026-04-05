# Scope

## In scope for the paper

**Variant types:**

- SNV, indel, and all SV types: fusions, translocations, large deletions, inversions, duplications, and tandem duplications.

**Evaluation modes:**

- Tumour-informed (monitoring) mode is the primary evaluation focus. This is where kam achieves strongest precision and is most clinically relevant.
- Discovery mode results are included but not the headline.

**Chemistry:**

- Configurable chemistry support via config.toml. Not hardcoded to Twist.
- Benchmarking covers public datasets beyond Twist to demonstrate generalisation.

**Comparisons:**

- Speed comparison against alignment-based methods.
- Accuracy comparison against thesis alignment-based results (RaSCALL) as the primary baseline.

## In scope for the tool (beyond paper)

- Production-ready CLI with `kam run --config config.toml` as the main entry point.
- Full configurability: chemistry, targets, thresholds, output options, filtering.
- CLI flags for overrides on top of config.toml.

## Out of scope

- Claiming superiority over alignment-based methods in all scenarios.
- De novo discovery as a primary selling point. It is secondary.
- Nextflow integration. Deferred past the paper.
- GUI or web interface.
- Multi-sample joint calling.
