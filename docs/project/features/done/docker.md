# Docker Image

## Status

In Progress.

## Summary

Multi-stage Dockerfile producing a slim runtime image with the kam binary and
ML models bundled. No reference genome needed. A `docker-compose.yml` file
provides easy invocation with bind-mounted data directories.

## Problem

No containerised deployment existed. Collaborators and reviewers who wanted to
run kam needed to install the Rust toolchain, clone the repository, and build
from source. This created friction for adoption and reproducibility,
particularly for reviewers evaluating benchmark results.

## Solution

A multi-stage Dockerfile at the repository root:

1. **Builder stage** (`rust:1.80-slim`): installs build dependencies, copies
   the workspace, and runs `cargo build --release`.
2. **Runtime stage** (`debian:bookworm-slim`): copies only the compiled `kam`
   binary from the builder. Installs minimal runtime dependencies (libc). ML
   models are already compiled into the binary via `include_bytes!`, so no
   extra model files are needed.

A `docker-compose.yml` at the repository root defines a `kam` service with
bind mounts for input data and output results.

The image does not include a reference genome. kam is alignment-free, so no
genome is required at runtime.

### Image size

The runtime image is expected to be under 100 MB. The kam binary is
approximately 30 MB (including bundled ML models). The base
`debian:bookworm-slim` image is approximately 60 MB.

## Tests

- `docker build -t kam:latest .` succeeds without errors
- `docker run --rm kam:latest --help` prints the usage message
- `docker-compose run kam --help` prints the usage message
- End-to-end: mount test FASTQs and targets, verify `variants.tsv` is produced

## Relevant Files

- `Dockerfile` — multi-stage build definition
- `docker-compose.yml` — compose service definition
- `.dockerignore` — excludes target/, bigdata/, .git/ from build context
