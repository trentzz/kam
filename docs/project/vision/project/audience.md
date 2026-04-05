# Audience

## Primary: Bioinformatics Researchers

Researchers interested in alignment-free methods and k-mer-based variant detection. They read the paper to understand the approach, its tradeoffs, and where it holds up.

They care about:

- Does the method work? What are the failure modes?
- How does it compare to alignment-based methods?
- What does molecule-level tracking buy you?

## Secondary: Clinical Bioinformaticians

Teams running ctDNA monitoring panels, particularly with Twist Biosciences duplex UMI chemistry. They want a fast, reliable pipeline they can drop into an existing workflow.

They care about:

- Sensitivity and precision on real samples.
- Reproducibility and deterministic output.
- Simple configuration without code changes.

## Tertiary: Tool Users

Bioinformaticians who want a fast, simple alternative to the HUMID + Jellyfish + km stack. They want a single binary with a config file.

They care about:

- Speed relative to the current pipeline.
- Correct output in the common case.
- Clear documentation for chemistry configuration.

## Paper readers vs tool users

Paper readers ask: does it work, how does it compare, and what are the tradeoffs?

Tool users ask: can I configure it for my chemistry, is it fast, and is the output trustworthy?

Both questions matter. The paper makes the case. The tool delivers it.
