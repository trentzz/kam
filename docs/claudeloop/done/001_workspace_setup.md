# Task 001: Rust Workspace Setup

## What to implement

Initialize the Rust workspace with all crate stubs. Each crate should compile with an empty `lib.rs` (or `main.rs` for the binary crate). Set up the workspace `Cargo.toml` with all member crates and shared dependency versions.

## Steps

1. Create root `Cargo.toml` as workspace with members: `kam-core`, `kam-assemble`, `kam-index`, `kam-pathfind`, `kam-call`, `kam`
2. Create each crate directory with `Cargo.toml` and minimal `src/lib.rs` (or `src/main.rs` for `kam`)
3. Add shared dependencies to workspace `Cargo.toml` under `[workspace.dependencies]`:
   - `thiserror = "2"`
   - `serde = { version = "1", features = ["derive"] }`
   - `log = "0.4"`
4. Each crate's `Cargo.toml` should reference workspace dependencies where applicable
5. `kam-assemble`, `kam-index`, `kam-pathfind`, `kam-call` all depend on `kam-core`
6. `kam` (binary) depends on all other crates

## Definition of done

- `cargo build` succeeds for the entire workspace
- `cargo test` passes (no tests yet, but no errors)
- `cargo clippy -- -D warnings` passes
- Directory structure matches `docs/planning/rust_workspace_architecture.md`
