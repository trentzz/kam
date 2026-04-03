# Paper Structure

## Narrative Arc

Alignment-based variant detection for ctDNA monitoring is slow and introduces reference bias. This paper asks whether alignment-free methods can match performance. The answer is yes, in monitoring mode, and the evidence is presented section by section.

---

## Sections

### 1. Introduction

The problem: alignment-based variant detection for ctDNA monitoring is slow and introduces reference bias. The question: can alignment-free methods match performance in tumour-informed monitoring?

State the core claim. Outline the paper.

### 2. Background

- Duplex UMI sequencing chemistry: Twist read structure, UMI design, molecule assembly.
- Current pipeline: HUMID, Jellyfish, km. What each stage does.
- Information loss at each stage. Why molecule-level tracking matters.

### 3. Related Work

- Alignment-based UMI tools: UMI-tools, fgbio, UMICollapse.
- k-mer tools without molecule awareness: Jellyfish, km.
- Position kam in the landscape: alignment-free, molecule-aware, chemistry-configurable.

### 4. Method

- kam architecture: five crates (kam-core, kam-assemble, kam-index, kam-pathfind, kam-call).
- Molecule assembly from raw FASTQ.
- k-mer indexing with molecule provenance.
- de Bruijn graph construction and path walking.
- Statistical variant calling.
- `config.toml` for chemistry configurability.
- Tumour-informed mode versus discovery mode.

### 5. Experimental Setup

- Datasets: Twist titration (5/15/30 ng, 0-2% VAF) and public datasets.
- Alignment-based baseline: RaSCALL results from the same samples.
- Metrics: sensitivity, precision, F1, runtime.
- Per-variant detailed reporting.

### 6. Results

- **SNV/indel**: tumour-informed sensitivity versus alignment-based across VAF levels.
- **SV**: detection across types (fusions, translocations, large deletions, inversions, duplications).
- **Speed**: runtime comparison on matched datasets.
- **Cross-chemistry**: results on non-Twist data.
- **Per-variant analysis**: where kam wins, where it loses, and why.

### 7. Discussion

- Tradeoffs between alignment-free and alignment-based approaches.
- When alignment-free is the right choice.
- Limitations: discovery mode sensitivity gaps, very low VAF.
- Future directions: production hardening, more chemistries, de novo discovery improvements.

### 8. Conclusion

Alignment-free variant detection is viable. In monitoring mode, kam matches or beats alignment-based methods. The approach generalises across chemistries and variant types.
