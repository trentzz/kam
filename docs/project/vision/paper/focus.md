# Paper Focus

## Core Claim

In tumour-informed monitoring mode, kam achieves comparable sensitivity to alignment-based methods with perfect precision, faster runtime, and no reference bias.

## Secondary Claims

- The approach generalises across UMI chemistries, not only Twist.
- SV detection (fusions, translocations, large deletions, inversions, duplications) is viable in an alignment-free framework.

## Central Story

When you know what you are looking for (monitoring mode), alignment-free is as good or better, faster, and simpler. It works across chemistries and variant types.

## What to Emphasise

- Tumour-informed monitoring results: precision 1.0, zero false positives.
- Speed advantage over alignment-based pipelines.
- Molecule-level evidence tracking as a differentiator from read-level methods.
- Configurability across UMI chemistries via `config.toml`.
- Broad SV type support in an alignment-free setting.

## What to Downplay

- Discovery mode sensitivity gaps: present these honestly, but they are not the headline.
- Limitations at very low VAF in discovery mode.

## What Not to Claim

- That alignment-free is universally better than alignment-based.
- That discovery mode is production-ready for clinical use.
