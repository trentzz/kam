# Task 027: Wire Up `kam call` Subcommand

## Location
`kam/src/commands/call.rs`

## What to implement

Connect the `kam call` CLI subcommand: read scored paths → call variants → write output in requested format(s) + QC.

## Interface

```rust
use crate::cli::CallArgs;

pub fn run_call(args: CallArgs) -> Result<(), Box<dyn std::error::Error>>;
```

## Logic

1. Read scored paths from bincode
2. For each target's scored paths:
   a. Identify reference path and variant paths
   b. Call variants using `call_variant`
3. Parse --output-format (comma-separated list)
4. Write variants to output file(s) using the requested format(s)
5. Build and write `CallQc`

## Tests

Integration test with temp files.

## Definition of done

- `cargo test -p kam` passes
- `cargo clippy -p kam -- -D warnings` passes
