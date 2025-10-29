# kam (Kmer Analysis Modules)

This is a collection of tools written as part of a pipeline for Kmer based variant detection.

## Docker

### From dockerhub

```bash
docker pull trentzz/kam
```

### From repo

Start by cloning this repo, then.

```bash
cd kam
docker build -t kam .
docker run -it --rm kam bash
```

You should now have access to the tools below, as well as jellyfish and km.

## Tools

| Tool                                                | Description                                                                                                                          |
| --------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------ |
| [multiseqex](https://github.com/trentzz/multiseqex) | **MULTI SEQuence EXtraction**. Batch process extracting sequences from a reference file. Similar to `samtools faidx` but multi-core. |
| [kmtools](https://github.com/trentzz/kmtools)       | Extension tools for km including multithreading, filtering, stats, and plotting.                                                     |
| [vcf2pandas](https://github.com/trentzz/vcf2pandas) | `vcf2pandas` is a python package to convert vcf files to `pandas` dataframes.                                                        |
| [vcf2xlsx](https://github.com/trentzz/vcf2xlsx)     | `vcf2xlsx` is a python cli package to convert vcf files to excel files.                                                              |
| [refolder](https://github.com/trentzz/refolder)     | CLI that moves files matching a glob pattern into an equal number of subfolders.                                                     |

## Workflow

You have a few options for using all the tools here and running a workflow:

1. Simple bash script: [example-workflow.sh](./example-workflow.sh)
2. Nextflow script: [example-workflow.nf](./example-workflow.nf)
3. Manually running each tool individually. This is useful when using a custom workflow, debugging, rerunning certain sections, etc.

### Bash Script

### Nextflow

### Manually
