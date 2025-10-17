# Base: Python 3.12 (Debian-based)
FROM python:3.12-slim

# Prevent interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    autoconf \
    automake \
    libtool \
    zlib1g-dev \
    curl \
    git \
    && rm -rf /var/lib/apt/lists/*

# Install Rust (latest stable)
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | bash -s -- -y
ENV PATH="/root/.cargo/bin:${PATH}"

# Install Jellyfish from source
RUN git clone https://github.com/gmarcais/Jellyfish.git /tmp/jellyfish && \
    cd /tmp/jellyfish && \
    autoreconf -i && \
    ./configure && \
    make && make install && \
    cd / && rm -rf /tmp/jellyfish

# Install pipx and set up its path
RUN pip install --no-cache-dir pipx && \
    pipx ensurepath
ENV PATH="/root/.local/bin:${PATH}"

# Install Python CLI tools with pipx
RUN pipx install km-walk && \
    pipx install vcf2xlsx

# Optional: verify installation
RUN jellyfish --version && \
    rustc --version && \
    km --help && \
    vcf2xlsx --help || true

# Default shell
CMD ["/bin/bash"]
