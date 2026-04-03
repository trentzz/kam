# Feature Priorities

## High priority — needed for paper

1. **config.toml support**: Full run configuration in one file. Chemistry parameters, targets, thresholds, output options, and filtering. CLI flags for overrides. `kam run --config config.toml` is the main entry point.

2. **SV type expansion**: Detection for fusions, translocations, large deletions, inversions, duplications, and tandem duplications in kam-pathfind.

3. **Chemistry configurability**: UMI length, position, duplex vs simplex, skip bases, all via config.toml. No hardcoded Twist assumptions.

4. **Public dataset benchmarking**: Broad search for public datasets across diverse chemistries and variant types. Demonstrates generalisation.

## Medium priority — important but not blocking

5. **Alignment-based comparison pipeline**: Systematic comparison against thesis RaSCALL results. Documents the baseline clearly.

6. **Per-variant detailed results**: Maintain current per-sample detail level for all new benchmarks.

## Lower priority — nice to have for paper

7. **Discovery mode improvements**: Improve sensitivity and reporting for de novo discovery.

8. **Additional output formats**: VCF improvements, MAF support, or other common formats.
