# Docker Usage

## Overview

The kam Docker image provides a self-contained environment for running the full variant detection
pipeline. The image contains:

- The `kam` binary with all built-in ML models bundled
- Minimal runtime dependencies (glibc)

The image does not contain a reference genome. kam is alignment-free, so no genome is needed at
any stage of the pipeline. All inputs (FASTQs, targets, junction sequences) are mounted at
runtime via bind mounts.

---

## Quick start

Pull the image from Docker Hub:

```bash
docker pull trentzz/kam:latest
```

Then run:

```bash
docker run --rm \
  -v $(pwd)/data:/data \
  -v $(pwd)/results:/results \
  trentzz/kam:latest run \
  --r1 /data/R1.fastq.gz \
  --r2 /data/R2.fastq.gz \
  --targets /data/panel.fa \
  --output-dir /results/
```

Results are written to `./results/` on the host.

---

## Building from source

Clone the repository and build the image:

```bash
git clone https://github.com/trentzz/kam
cd kam
docker build -t kam:latest .
```

The build uses a multi-stage Dockerfile:

1. **Builder stage** (`rust:1.94-slim`): compiles the workspace in release mode.
2. **Runtime stage** (`debian:trixie-slim`): copies only the compiled binary. No Rust
   toolchain or source code is present in the final image.

Build time is approximately 5-10 minutes on a modern machine. The final image is under 100 MB.

---

## Using docker-compose

The repository includes a `docker-compose.yml` at the root. It defines a `kam` service with
bind mounts for common directories.

```bash
# Run the full pipeline
docker compose run --rm kam run \
  --r1 /data/R1.fastq.gz \
  --r2 /data/R2.fastq.gz \
  --targets /data/panel.fa \
  --output-dir /results/

# List available models
docker compose run --rm kam models list
```

The compose file mounts `./data` to `/data` and `./results` to `/results` inside the container.
Place your input files in `./data` and find outputs in `./results`.

To use different host directories, override the mount paths:

```bash
docker compose run --rm \
  -v /path/to/my/fastqs:/data \
  -v /path/to/my/output:/results \
  kam run \
  --r1 /data/sample_R1.fq.gz \
  --r2 /data/sample_R2.fq.gz \
  --targets /data/panel.fa \
  --output-dir /results/
```

---

## Passing configuration files

Mount a directory containing your TOML config file and reference it with `--config`:

```bash
docker run --rm \
  -v $(pwd)/data:/data \
  -v $(pwd)/results:/results \
  -v $(pwd)/config:/config \
  trentzz/kam:latest run \
  --config /config/kam.toml \
  --output-dir /results/
```

The config file paths must reference the container paths, not the host paths. If the config
file contains relative paths for `r1`, `r2`, or `targets`, they resolve relative to the
working directory inside the container. Use absolute container paths for clarity:

```toml
# kam.toml — paths as seen inside the container
[input]
r1 = "/data/sample_R1.fq.gz"
r2 = "/data/sample_R2.fq.gz"
targets = "/data/panel.fa"

[output]
output_dir = "/results/"
output_format = "tsv,vcf"
```

---

## Junction sequences from BAM

To monitor a fusion or SV observed in a BAM file, mount the junction FASTA alongside your
other inputs:

```bash
docker run --rm \
  -v $(pwd)/data:/data \
  -v $(pwd)/results:/results \
  trentzz/kam:latest run \
  --r1 /data/plasma_R1.fastq.gz \
  --r2 /data/plasma_R2.fastq.gz \
  --targets /data/panel.fa \
  --junction-sequences /data/junction_sequences.fa \
  --output-dir /results/ \
  --output-format tsv,vcf
```

See the [patient SV monitoring guide](guides/patient-sv-monitoring.md) for how to create the
junction FASTA from an IGV observation.

---

## SV and fusion detection

Mount additional input files for SV junctions or fusion targets:

```bash
docker run --rm \
  -v $(pwd)/data:/data \
  -v $(pwd)/results:/results \
  trentzz/kam:latest run \
  --r1 /data/R1.fastq.gz \
  --r2 /data/R2.fastq.gz \
  --targets /data/panel.fa \
  --sv-junctions /data/sv_junctions.fa \
  --fusion-targets /data/fusions.fa \
  --output-dir /results/ \
  --output-format tsv,vcf
```

---

## Output files

All outputs go to the directory specified by `--output-dir`, which maps to the mounted
`/results` volume. On the host, find:

| File | Contents |
|---|---|
| `results/variants.tsv` | All variant calls (PASS and filtered) |
| `results/variants.vcf` | VCF output (if `--output-format` includes `vcf`) |
| `results/assembly_qc.json` | Molecule assembly statistics |
| `results/index_qc.json` | K-mer indexing statistics |
| `results/pathfind_qc.json` | Path finding statistics |
| `results/call_qc.json` | Variant calling statistics |

---

## Tumour-informed monitoring

```bash
docker run --rm \
  -v $(pwd)/data:/data \
  -v $(pwd)/results:/results \
  trentzz/kam:latest run \
  --r1 /data/plasma_R1.fastq.gz \
  --r2 /data/plasma_R2.fastq.gz \
  --targets /data/panel.fa \
  --target-variants /data/tumour_biopsy.vcf \
  --output-dir /results/ \
  --max-vaf 0.35 \
  --output-format tsv,vcf
```

---

## ML model selection

Built-in models are bundled inside the binary. Select one by name:

```bash
docker run --rm \
  -v $(pwd)/data:/data \
  -v $(pwd)/results:/results \
  trentzz/kam:latest run \
  --r1 /data/R1.fastq.gz \
  --r2 /data/R2.fastq.gz \
  --targets /data/panel.fa \
  --ml-model twist-duplex-v2 \
  --output-dir /results/
```

To use a custom ONNX model, mount the model file:

```bash
docker run --rm \
  -v $(pwd)/data:/data \
  -v $(pwd)/results:/results \
  -v $(pwd)/models:/models \
  trentzz/kam:latest run \
  --r1 /data/R1.fastq.gz \
  --r2 /data/R2.fastq.gz \
  --targets /data/panel.fa \
  --custom-ml-model /models/my_model.onnx \
  --output-dir /results/
```

The companion `.json` metadata file must be in the same directory as the ONNX file.

---

## Resource requirements

| Resource | Recommendation |
|---|---|
| CPU | 1 core minimum. kam currently runs single-threaded. |
| RAM | 2 GB for standard panels (< 500 targets) at 2M reads. 4 GB for large panels or high-depth samples. |
| Disk | 500 MB for the image. Output files are typically < 10 MB. |

RAM usage scales with panel size (number of targets) and sequencing depth (number of unique
molecules). The indexing stage holds all k-mers in memory. For very large panels (> 1000
targets), allocate 4-8 GB.

Set memory limits in Docker if running multiple containers:

```bash
docker run --rm --memory=4g \
  -v $(pwd)/data:/data \
  -v $(pwd)/results:/results \
  trentzz/kam:latest run \
  --r1 /data/R1.fastq.gz \
  --r2 /data/R2.fastq.gz \
  --targets /data/panel.fa \
  --output-dir /results/
```

---

## Batch processing

Process multiple samples by running one container per sample:

```bash
for sample in sample1 sample2 sample3; do
  docker run --rm \
    -v $(pwd)/data:/data \
    -v $(pwd)/results:/results \
    trentzz/kam:latest run \
    --r1 "/data/${sample}_R1.fastq.gz" \
    --r2 "/data/${sample}_R2.fastq.gz" \
    --targets /data/panel.fa \
    --output-dir "/results/${sample}/" \
    --output-format tsv,vcf
done
```

Each run creates its output directory inside the mounted `/results` volume.

---

## Troubleshooting

### apt-get fails during build

If `docker build` fails with "Unable to locate package", Docker's default bridge networking may not be resolving DNS correctly. Build with host networking:

```bash
docker build --network host -t trentzz/kam:latest .
```

This is common when a VPN is active.
