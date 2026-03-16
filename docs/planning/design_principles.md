# kam Design Principles

## One Word: Transparent

Every decision kam makes — every read dropped, family merged, variant called or filtered — should be explainable, traceable, and configurable. The user should never be left wondering "why did it do that?"

## Core Principles

### 1. Granular User Control

Every threshold, heuristic, and behavioural choice should be configurable. kam should never force a tradeoff that only the user can evaluate.

- Every filter has a parameter. No magic numbers buried in source code.
- Sensible defaults for common use cases (Twist duplex ctDNA panels), but defaults are always overridable.
- If someone wants to try a different approach with this tool, they should be able to. The tool serves the user's hypothesis, not the developer's assumptions.

**Examples:**
- `--min-umi-quality 20` (default Q20, but user can set Q10 for exploratory work or Q30 for high-stringency)
- `--fingerprint-max-diff 4` (default 4, but tighten to 2 for high-quality data or loosen to 6 for FFPE)
- `--min-family-size 3` / `--min-duplex-reads 1` / etc.
- Chemistry presets (`--chemistry twist-umi-duplex`) set many defaults at once, but each can be individually overridden

### 2. Extreme Logging and Observability

The user should be able to understand exactly what happened during a run, at whatever level of detail they need.

- **Always-on counters:** summary statistics are always collected (negligible cost). Every run produces a `summary.json` with counts of what happened.
- **Per-category log files:** each type of event (dropped reads, UMI clustering decisions, filtered variants, etc.) has its own log file, independently togglable.
- **Zero cost when off:** disabled logs compile to nothing. Production runs are fast. Debug runs are detailed.
- **`--log all`** for full transparency, **`--log none`** for maximum speed.

See `docs/planning/logging_architecture.md` for the full logging design.

### 3. Explain Why, Not Just What

When kam drops a read, filters a variant, or merges two UMI groups, the log should say *why*.

- "Read READ_042 dropped: r1 UMI base at position 3 has quality Q8, below threshold Q20"
- "Variant at TP53_exon7 filtered: STRAND_BIAS (p=0.003, all 4 supporting molecules are forward-strand only)"
- "UMI group ACGTA+TGCAT split into 2 families: endpoint fingerprints differ by 38 bits (threshold: 4)"

This is not just for debugging — it's for scientific understanding. When a clinician asks "why wasn't this variant called?" or "why was this variant called?", the answer should be in the logs.

### 4. No Black Boxes

- The statistical model is documented and the math is accessible. Users should understand what P(variant) means and what assumptions it relies on.
- QC metrics have clear definitions. "duplex_fraction: 0.798" means exactly "79.8% of molecule families have at least one read on each strand."
- Filter labels (PASS, STRAND_BIAS, LOW_CONFIDENCE, LOW_DUPLEX, COLLISION_RISK) have documented criteria that reference specific thresholds the user can change.

### 5. Configurable Strictness Spectrum

Different use cases need different tradeoffs. kam should support the full spectrum:

| Use Case | Strictness | Example Config |
|----------|-----------|----------------|
| Clinical MRD monitoring | Maximum | High family size, duplex required, strict UMI quality |
| Research/exploratory | Moderate | Lower thresholds, simplex accepted, report everything |
| Method development | Minimal filtering | Accept all reads, log everything, no filters |
| Benchmarking | Full pipeline | Standard defaults with sensitivity tracking |

Rather than separate "modes," this is achieved through the parameter space. Chemistry presets provide good starting points.

### 6. Reproducibility

- All output is deterministic for given inputs and parameters
- Seeded RNG everywhere (for any stochastic operations)
- Parameters are recorded in output headers and QC JSON
- Exact command line and version captured in every output file
- Nextflow resume/caching requires this — it's not optional

### 7. Fail Loud, Not Silent

- If something is wrong, error immediately with a clear message. Never silently produce bad results.
- QC checks between pipeline stages: each stage validates the previous stage's output before proceeding.
- Warnings for suspicious metrics (e.g., duplex fraction below 50%, collision rate above 5%)
- Non-zero exit codes for any failure, with structured error output that Nextflow can act on.

## How These Principles Apply to Implementation

### For every parameter added:
1. Document what it controls and why you'd change it
2. Set a sensible default with a comment explaining the reasoning
3. Include it in the CLI help text with a concrete example
4. Record the value used in QC JSON output
5. Log when it causes an action (read dropped, variant filtered, etc.)

### For every heuristic or decision point:
1. Make it a named, configurable threshold (not a hardcoded constant)
2. Log the decision and the values that led to it (when logging enabled)
3. Include the threshold value in filter labels so the user can see what triggered it

### For every output:
1. Include a header or metadata section with the parameters used
2. Machine-readable format (TSV, JSON) alongside human-readable summaries
3. Enough detail to reproduce the run and understand every line of output
