# Base: Python 3.12 (Debian-based)
FROM python:3.12.11-slim-bookworm

# Prevent interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

# Install system dependencies (including build tools)
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    autoconf \
    automake \
    libtool \
    pkg-config \
    zlib1g-dev \
    curl \
    git \
    ca-certificates \
    openjdk-17-jdk \
    && rm -rf /var/lib/apt/lists/*

# Install Rust (latest stable)
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | bash -s -- -y
ENV PATH="/root/.cargo/bin:${PATH}"

# Install setuptools for Python 3.12 (needed for legacy distutils import)
RUN pip install --no-cache-dir setuptools wheel

# Download and install jellyfish
RUN curl -L https://github.com/gmarcais/Jellyfish/releases/download/v2.2.6/jellyfish-2.2.6.tar.gz --output jellyfish-2.2.6.tar.gz
RUN tar zxvf jellyfish-2.2.6.tar.gz && \
    cd jellyfish-2.2.6 && \
    ./configure PYTHON=$(which python3.12) --prefix=$VIRTUAL_ENV --enable-python-binding && \
    make -j 4 && \
    make install 

# Clean up
RUN rm -rf jellyfish-2.2.6 && \
    rm -rf jellyfish-2.2.6.tar.gz

# Install pipx and ensure path
RUN pip install --no-cache-dir setuptools wheel pipx && \
    pipx ensurepath
ENV PATH="/root/.local/bin:${PATH}"

# Install Python CLI tools with pipx
RUN pipx install km-walk && \
    pipx install vcf2xlsx && \
    pip install vcf2pandas && \
    pipx install git+https://github.com/trentzz/kmtools

# Install Rust-based tools via cargo
RUN cargo install --git https://github.com/trentzz/multiseqex && \
    cargo install refolder

# Install Nextflow
RUN curl -s https://get.nextflow.io | bash && \
    mv nextflow /usr/local/bin/ && \
    chmod +x /usr/local/bin/nextflow

# Verify installations (non-fatal)
RUN jellyfish --version && \
    rustc --version && \
    km --help && \
    multiseqex --help && \
    refolder --help && \
    kmtools --help && \
    vcf2xlsx --help || true

# Default shell
CMD ["/bin/bash"]
