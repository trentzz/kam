# Testing Review

## Test Statistics

- **Total tests**: 293 (all passing)
- **Unit tests**: ~250 (in `#[cfg(test)] mod tests` blocks)
- **Doc tests**: ~43 (in `///` examples)
- **Integration tests**: 2 (in `kam/src/commands/run.rs`)
- **Clippy**: Zero warnings with `-D warnings`
- **Test execution time**: ~2 seconds total

## Strengths

### Every Function Has a Test

This is not hyperbole — I checked. Every public function in every crate has at least one corresponding test. Most have 3-5 tests covering normal operation, edge cases, and error conditions. The test-per-function discipline is exceptional.

### Doc Tests That Actually Run

The doc examples are compiled and executed as tests. They catch drift between documentation and implementation. With 43 doc tests, this is a meaningful safety net.

### Test Naming Convention

Tests follow a consistent descriptive naming pattern:
- `valid_read_pair_long_template`
- `strand_bias_balanced_is_nonsignificant`
- `cluster_no_transitive_chaining`

These read as specifications. You can understand what the function should do just by reading the test names.

### Good Edge Case Coverage

The tests include meaningful edge cases:
- Empty inputs (empty FASTQ, empty molecule list, empty k-mer set)
- Boundary conditions (exactly at threshold, exactly at minimum length)
- Degenerate cases (all bases masked, single-read family, self-loop k-mer)
- Error conditions (mismatched FASTQ counts, short reads, bad magic bytes)

### Determinism Testing

Several tests verify deterministic output:
- `fingerprint_deterministic` — same input produces same fingerprint across 10 calls
- `canonical_pair_hash_consistent` — hash set membership is stable
- `assembly_stats_counts_correct` — specific count expectations

## Weaknesses

### 1. No Real Data Tests

Every test uses synthetic data constructed inline. There are no tests with:
- Real FASTQ files (even small ones)
- Real Twist UMI chemistry output
- Known variant positive controls
- Known negative controls (wildtype samples)

This means the code has never been validated against the actual data format it's designed to process. Real FASTQ has edge cases that synthetic data doesn't: adapter contamination, polyG tails, quality score distributions, read name formats.

**Recommendation**: Add a `tests/fixtures/` directory with a small real or realistic FASTQ pair (100 read pairs) and at least one integration test that runs the full pipeline on it.

### 2. No Property-Based Testing

For a pipeline with mathematical invariants, property-based testing (using `proptest` or `quickcheck`) would catch edge cases that hand-written tests miss:

- `canonical(kmer) == canonical(reverse_complement(kmer))` for all k-mers
- `umi_pair_hamming_distance(a, b) == umi_pair_hamming_distance(b, a)` for all pairs
- `encode_kmer(decode_kmer(x, k), k) == x` for all valid `x, k`
- `estimate_vaf(k, m).1 <= estimate_vaf(k, m).0 <= estimate_vaf(k, m).2` for all `k <= m`

These invariants are true by construction, but property tests would verify them across millions of random inputs.

### 3. No Regression Tests

There's no mechanism to catch regressions in numerical output. If the VAF estimator changes slightly due to a floating-point optimization, no test would catch it. Consider adding golden-file tests: run the pipeline on fixed input, compare output byte-for-byte against a checked-in expected output.

### 4. Consensus Tests Don't Test the Hard Cases

The consensus calling tests verify:
- All reads agree → correct base called
- Majority wins over minority
- Low-quality bases excluded
- Duplex agreement multiplies error probs

But they don't test:
- What happens with 3 different bases at one position (A, C, T)?
- What happens when all reads are below quality threshold?
- What happens with extremely unbalanced families (99 fwd, 1 rev)?
- What's the consensus of a single base when it's N?

### 5. Walk Tests Use Artificial Graph IDs

The path walking tests construct graphs with integer k-mer IDs (0, 1, 2, ...) rather than real encoded k-mers:

```rust
let mut g = DeBruijnGraph::new(k);
g.add_edge(0, 1);
g.add_edge(1, 2);
```

This tests the graph traversal logic but doesn't test the k-mer overlap logic that `from_index()` uses. The `from_index()` tests in `graph.rs` do use real k-mers, which is good, but the walk tests are decoupled from the encoding.

### 6. Fisher's Exact Test Not Validated Against R

The strand bias test implements Fisher's exact test from scratch. The implementation looks correct, but it's not validated against a reference implementation (R's `fisher.test()`, scipy's `fisher_exact()`). For a statistical test that gates clinical variant calls, this should be cross-validated.

### 7. No Fuzzing

For a tool that parses untrusted input (FASTQ files from sequencers), there should be fuzz tests for:
- The FASTQ parser (malformed records, truncated files)
- The bincode deserializer (corrupted intermediate files)
- The FASTA reader (unusual line lengths, non-ASCII characters)

### 8. Test Isolation

Several tests use `std::env::temp_dir()` with PID-based names:

```rust
fn tmp_path(name: &str) -> std::path::PathBuf {
    let mut p = std::env::temp_dir();
    p.push(format!("kam_core_serialize_test_{name}_{}", std::process::id()));
    p
}
```

The `io.rs` tests use a custom `TempDir` RAII struct instead of the `tempfile` crate that's already a dependency (used in `run.rs` tests). This inconsistency should be unified — just use `tempfile::tempdir()` everywhere.

### 9. No Benchmarks

Despite `docs/benchmarking/` being listed in CLAUDE.md, the directory is empty and no `cargo bench` targets exist. For a tool that claims to replace Jellyfish (a highly optimized C tool), performance testing is essential.

## Test Quality by Crate

| Crate | Tests | Quality | Gap |
|---|---|---|---|
| kam-core | 21 | Good | Could use property tests for molecule types |
| kam-assemble | 62 | Very Good | Missing real-data tests, edge cases in consensus |
| kam-index | 34 | Good | No large-scale indexing tests |
| kam-pathfind | 40 | Good | Walk tests use artificial graph IDs |
| kam-call | 32 | Good | Fisher test needs cross-validation |
| kam (CLI) | 18 | Good | Only 2 integration tests |

## Recommendations

1. **Add fixture-based integration tests** with realistic (even if small) FASTQ data
2. **Add `proptest` for mathematical invariants** (encoding roundtrips, symmetries)
3. **Cross-validate Fisher's exact test** against R/scipy
4. **Add golden-file regression tests** for the full pipeline
5. **Unify temp file handling** on `tempfile` crate
6. **Add `criterion` benchmarks** for hot paths (k-mer encoding, consensus, graph construction)
