# Dockerfile for kam — alignment-free variant detection for duplex UMI sequencing.
#
# Multi-stage build:
#   Stage 1 (builder): compiles the full Rust workspace with ML support enabled.
#   Stage 2 (runtime): minimal Debian image containing only the binary and
#                      the ONNX Runtime shared library used by ML inference.
#
# kam is alignment-free and works directly from FASTQ + target FASTA.
# No reference genome is bundled or required.
#
# Build:
#   docker build -t kam:0.3.0 .
#
# Run:
#   docker run --rm \
#     -v $(pwd)/data:/data \
#     -v $(pwd)/results:/results \
#     kam:0.3.0 run \
#       --r1 /data/R1.fastq.gz \
#       --r2 /data/R2.fastq.gz \
#       --targets /data/panel.fa \
#       --output-dir /results/

# ---------------------------------------------------------------------------
# Stage 1: Build
# Uses the official Rust image to compile the workspace with ML support.
# ---------------------------------------------------------------------------
FROM rust:1.94-slim AS builder

# Install build-time dependencies needed by openssl and other sys crates.
RUN apt-get update && apt-get install -y --no-install-recommends \
    pkg-config \
    libssl-dev \
    g++ \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /build

# Copy manifests first so that cargo can cache the dependency layer.
# When only source files change, this layer is reused.
COPY Cargo.toml Cargo.lock ./
COPY kam-core/Cargo.toml         kam-core/Cargo.toml
COPY kam-assemble/Cargo.toml     kam-assemble/Cargo.toml
COPY kam-index/Cargo.toml        kam-index/Cargo.toml
COPY kam-pathfind/Cargo.toml     kam-pathfind/Cargo.toml
COPY kam-call/Cargo.toml         kam-call/Cargo.toml
COPY kam-ml/Cargo.toml           kam-ml/Cargo.toml
COPY kam/Cargo.toml              kam/Cargo.toml

# Create stub lib/main files for each crate so cargo can resolve and
# compile dependencies without the real source.
RUN mkdir -p kam-core/src && echo "pub fn _stub() {}" > kam-core/src/lib.rs && \
    mkdir -p kam-assemble/src && echo "pub fn _stub() {}" > kam-assemble/src/lib.rs && \
    mkdir -p kam-index/src && echo "pub fn _stub() {}" > kam-index/src/lib.rs && \
    mkdir -p kam-pathfind/src && echo "pub fn _stub() {}" > kam-pathfind/src/lib.rs && \
    mkdir -p kam-call/src && echo "pub fn _stub() {}" > kam-call/src/lib.rs && \
    mkdir -p kam-ml/src && echo "pub fn _stub() {}" > kam-ml/src/lib.rs && \
    mkdir -p kam/src && echo "fn main() {}" > kam/src/main.rs

# Compile dependencies only (source is stubs — changes to real source
# will not invalidate this layer).
RUN cargo build --release --features ml 2>/dev/null || true

# Now copy the real source and rebuild.
COPY kam-core/src     kam-core/src
COPY kam-assemble/src kam-assemble/src
COPY kam-index/src    kam-index/src
COPY kam-pathfind/src kam-pathfind/src
COPY kam-call/src     kam-call/src
COPY kam-ml/src       kam-ml/src
COPY kam/src          kam/src
COPY kam/models       kam/models

# Touch source files so cargo knows they changed after the stub build.
RUN find . -name "*.rs" -exec touch {} +

RUN cargo build --release --features ml

# ---------------------------------------------------------------------------
# Stage 2: Runtime
# Minimal Debian image — no Rust toolchain, no build tools.
# ---------------------------------------------------------------------------
FROM debian:trixie-slim AS runtime

LABEL org.opencontainers.image.title="kam" \
      org.opencontainers.image.description="Alignment-free variant detection for duplex UMI sequencing" \
      org.opencontainers.image.version="0.3.0" \
      org.opencontainers.image.source="https://github.com/trentzz/kam"

# Install runtime libraries required by the binary and the ONNX Runtime.
# Note: if build fails here with "Unable to locate package", try building
# with --network host (e.g. docker build --network host -t kam .)
RUN apt-get update && apt-get install -y --no-install-recommends \
    libssl3 \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# Copy the compiled binary from the builder stage.
COPY --from=builder /build/target/release/kam /usr/local/bin/kam

# Copy the ONNX Runtime shared library downloaded by the ort crate at build
# time (via the download-binaries feature). The library lands in the build
# output directory; we place it in /usr/local/lib so the dynamic linker finds
# it without requiring LD_LIBRARY_PATH.
#
# The glob covers both the unversioned and versioned symlink names that ort
# may produce (e.g. libonnxruntime.so, libonnxruntime.so.1.18.0).
COPY --from=builder /build/target/release/build/ort-sys-*/out/libonnxruntime* \
    /usr/local/lib/

# Rebuild the dynamic linker cache so /usr/local/lib is searched.
RUN ldconfig

# Run as a non-root user for safety.
RUN useradd --no-create-home --shell /usr/sbin/nologin kam
USER kam

ENTRYPOINT ["kam"]
