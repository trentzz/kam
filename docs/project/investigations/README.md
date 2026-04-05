# investigations

Diagnostic investigation writeups. Each document records a problem encountered during development: the symptom, what was tested, the root cause, and the resolution. These are permanent records, not active tasks.

## Index

| File | Date | Area | Symptom | Root cause | Resolution |
|---|---|---|---|---|---|
| `anchor_missing_investigation.md` | 2026-03-22 | SNV/indel, pathfind | 18% of targets had missing start anchors in the graph | Exact k-mer anchor absent from raw index for many indel targets | Soft anchor fallback implemented |
| `sensitivity_investigation.md` | 2026-03-21 | SNV/indel, call | Sensitivity 51–61% at 2% VAF | Multiple causes: anchor gaps, family size filtering, low-VAF dropout | Documented; guided subsequent improvements |
| `fn_investigation_2026-03-20.md` | 2026-03-20 | SNV/indel | FN analysis at 30ng 2% VAF | Low duplex fraction (6.6%) reduces evidence quality | Identified duplex-specific scoring gaps |
| `fp_investigation/` | 2026-03-21 | SNV/indel, call | False positive calls under various conditions | Multi-round investigation; resolved through filter tuning | Four rounds documented |
| `depth_saturation_investigation.md` | — | SNV/indel | Sensitivity plateaus above certain read depths | Coverage saturation above molecule family capacity | Documented read depth sweet spot |
| `sv_benchmark_investigation.md` | — | SV | SV benchmark setup and results analysis | — | Baseline established |
| `sv_threshold_investigation.md` | — | SV, call | Optimal threshold selection for SV calling | — | Per-type thresholds determined |
| `sv_large_detection_investigation.md` | — | SV | Detection of large SVs (>100 bp) | — | k-mer strategy for large SVs documented |
| `sv_sensitivity_report.md` | — | SV | SV sensitivity summary across types and VAFs | — | Sensitivity baseline recorded |
| `invdel_misclassification_investigation.md` | — | SV | INVDEL variants misclassified as other types | — | Classifier fix applied |
| `novins_misclassification_investigation.md` | — | SV | Novel insertion variants misclassified | — | Classifier fix applied |
| `inversion_detection.md` | — | SV | Inversion detection approach | — | Detection strategy documented |
| `dup_junction_detection.md` | — | SV | Tandem duplication junction detection | — | Detection strategy documented |
| `fusion_vaf_zero_investigation.md` | — | Fusion | Fusion calls with VAF = 0 | — | Root cause and fix documented |
| `target_length_investigation.md` | — | Pathfind | Effect of target sequence length on detection | — | Optimal target length range determined |
| `snvindel_benchmark_scoring_investigation.md` | — | SNV/indel | Scoring discrepancies in benchmark | — | Scoring logic corrected |
| `snvindel_results_v10_2m_reads_k31.md` | — | SNV/indel | Results for v10, 2M reads, k=31 | — | Full results table and figures |

## Writing a new investigation

Follow the template in `docs/project/devmanual/investigation-docs.md`.
