# devmanual

Developer reference for contributing to kam. Read the relevant file before starting work in that area.

## Contents

| File | When to read |
|---|---|
| `benchmarking.md` | Running benchmarks, adding datasets, scoring results |
| `task-tracking.md` | Picking up tasks, creating epics, committing work |
| `investigation-docs.md` | Writing up a diagnostic investigation |
| `release.md` | Preparing and publishing a release |
| `nextcloud.md` | Downloading and uploading large files via Nextcloud |
| `rust_workspace_architecture.md` | Crate layout, dependencies, build order |
| `core_data_model.md` | Molecule, ConsensusRead, MoleculeEvidence, and shared traits |
| `design_principles.md` | Transparency, granular control, logging, reproducibility |
| `logging_architecture.md` | Per-file togglable logs, drop log format, Nextflow interaction |
| `nextflow_integration.md` | HPC config, QC JSON, morning report |
| `development_workflow.md` | Claude Code loop, task format, safety boundaries |
| `current_pipeline_analysis.md` | Information loss analysis at each stage of the legacy pipeline |
| `fusion_detection_design.md` | Fusion detection design |
| `de_novo_discovery_design.md` | De novo discovery: 4 approaches, panel-aware design |
| `sv_detection_design.md` | SV detection design |
| `sv_detection_implementation.md` | SV detection implementation details |
| `output_format_specs.md` | Annotated FASTQ tag definitions, QC JSON schemas, TSV/VCF/JSON output |
| `hash_partition_umi_grouping.md` | Hash-partition UMI grouping algorithm |
| `ml_gradient_boosting.md` | Gradient boosting for variant re-scoring |
| `ml_integration_design.md` | ML integration design: ONNX runtime, feature extraction, bundled models |
