# Results Priorities

## Priority Order

1. **Tumour-informed sensitivity versus alignment-based at matched VAF levels.** This result makes or breaks the paper. It must show that kam is comparable.
2. **Precision comparison.** kam's zero false positive rate is a strong, clean result. Emphasise it.
3. **SV detection sensitivity across all supported types.** Demonstrates breadth.
4. **Runtime comparison.** Supports the speed advantage claim.
5. **Cross-chemistry generalisation.** Supports the scope of the approach.
6. **Discovery mode results.** Secondary. Present honestly without overstating.

---

## Metrics

For each experiment, report:

- **Sensitivity** (true positive rate) at each VAF level.
- **Precision** (positive predictive value).
- **F1 score.**
- **Runtime** in wall-clock seconds.
- **Per-variant detection concordance**: which specific variants does each method find?

---

## What Does Not Matter

- Memory usage, unless it becomes a practical bottleneck.
- Theoretical complexity analysis.
- Comparison against tools that do not handle UMI data.
