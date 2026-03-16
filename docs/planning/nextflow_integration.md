# Nextflow Integration

## Pipeline Structure

```groovy
// main.nf
nextflow.enable.dsl = 2

workflow {
    reads_ch = Channel.fromPath(params.samplesheet)
        .splitCsv(header: true)
        .map { row -> tuple(row.sample_id, file(row.fastq_r1), file(row.fastq_r2)) }

    molecules    = MOLECULE_ASSEMBLY(reads_ch)
    kmer_index   = KMER_INDEXING(molecules.output, params.targets)
    graph        = GRAPH_WALK(kmer_index.output, params.targets)
    calls        = VARIANT_CALLING(graph.output)
    MORNING_REPORT(calls.output.collect(), calls.qc.collect())
}
```

## Key Nextflow Features for kam

### Resume (Content-Addressed Caching)
`-resume` hashes each process execution based on inputs, script, and container. Identical inputs → cached output reused. **Critical requirement:** all kam stages must be deterministic for given inputs (no HashMap iteration order in output, no unseeded RNG, no timestamps).

### Error Strategy (Unattended HPC Running)
```groovy
process {
    errorStrategy = { task.exitStatus in [130, 137, 140, 143]
                      ? 'retry'    // OOM or walltime — retry with more resources
                      : 'finish'   // other errors — finish other jobs then fail
                    }
    maxRetries = 3
}
```

Exit codes 130, 137, 140, 143 = OOM kill / walltime exceeded on most HPC schedulers. Retry multiplies resources via `task.attempt`.

### Resource Profiles
```groovy
profiles {
    hpc {
        process.executor = 'slurm'
        withLabel: 'high_memory' { memory = { 32.GB * task.attempt }; cpus = 16 }
        withLabel: 'kmer_heavy'  { memory = { 16.GB * task.attempt }; cpus = 8 }
        withLabel: 'quick'       { memory = 4.GB; cpus = 4 }
    }
    local { process.executor = 'local'; process.cpus = 8; process.memory = 16.GB }
    test  { params.samplesheet = "test_data/samples_small.csv"; process.memory = 4.GB }
}
```

## Structured QC Output

Each kam stage outputs a QC JSON that the next stage validates:

```json
{
  "stage": "molecule_assembly",
  "sample_id": "PATIENT_001",
  "timestamp": "2026-03-14T02:17:43Z",
  "git_sha": "a3f8b2c1",
  "metrics": {
    "n_input_read_pairs": 4821044,
    "n_molecules": 612847,
    "duplex_fraction": 0.798,
    "mean_family_size": 7.2
  },
  "passed": true
}
```

Each stage runs `kam check-qc --qc-file <prev_stage_qc> --stage <name>` before proceeding. Non-zero exit = Nextflow sees failure and retries/fails cleanly.

## Morning Report

Final pipeline step aggregating all QC and variant calls:

```groovy
process MORNING_REPORT {
    publishDir "${params.outdir}/reports", mode: 'copy'
    input: path qc_files; path call_files; path trace
    script: """
    kam generate-report --qc-dir . --calls-dir . --nextflow-trace ${trace} \
        --output morning_report.txt --metrics sensitivity_metrics.tsv
    """
}
```

Uses Nextflow's own `trace.tsv` for wall times, peak memory, retry counts.

## Overnight Submission

```bash
nextflow run main.nf -profile hpc -resume --outdir "results/$(date +%Y%m%d)" \
    -params-file params/nightly.yaml
```

`-resume` always on. `-with-tower` (Seqera Platform free tier) for web dashboard monitoring without SSH.

## Sensitivity Tracking

Version-control `sensitivity_metrics.tsv` after every run:
```bash
cp results/*/reports/sensitivity_metrics.tsv metrics_history/$(date +%Y%m%d).tsv
git add metrics_history/ && git commit -m "nightly metrics: 0.1% VAF sensitivity XX%"
```

Git-tracked history of sensitivity at every VAF level, correlated with code changes.
