## Release Configuration

### Registry
- Primary: crates.io (cargo publish)
- Secondary: GitHub releases (tagged)

### Publish command
```
cargo publish -p kam-core
cargo publish -p kam-assemble
cargo publish -p kam-index
cargo publish -p kam-pathfind
cargo publish -p kam-call
cargo publish -p kam
```
Publish in dependency order. kam-core first, then crates that depend on it, then the binary crate last.

### Package contents
Rust source code only. Exclude:
- docs/ (benchmarks, paper, research, planning, vision, claudetracking, claudeguide)
- examples/ (config examples — these are for the repo, not the crate)
- Benchmark data and results
- Any non-source files

### Auth
crates.io token already configured.

### Pre-publish checks
```
cargo fmt -- --check
cargo clippy -- -D warnings
cargo test
```

### Post-publish steps
1. Create GitHub release from the tag
2. Attach the binary (if cross-compiled)

### Exclusions
Each crate's Cargo.toml should have:
```toml
[package]
exclude = ["docs/", "benchmarking/", "tests/fixtures/"]
```

The workspace root Cargo.toml should exclude docs/, examples/, and benchmark data from the published package.

### Notes
- Publish crates in dependency order to avoid "not found" errors
- Wait a few seconds between publishes for crates.io indexing
- The kam binary crate depends on all library crates
