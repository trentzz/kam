# Project Goals

## Purpose

Alignment-free variant detection for duplex UMI sequencing. kam replaces the HUMID + Jellyfish + km + kmtools pipeline with a unified Rust workspace that preserves molecule-level information throughout.

Target use case: detecting somatic variants at low VAF in liquid biopsy ctDNA panel samples.

## What success looks like

A published proof-of-concept paper showing the alignment-free approach is viable and comparable to alignment-based methods. A tool that works across chemistries and variant types, built with production quality.

## What this project IS

- A configurable, alignment-free variant detection pipeline.
- A proof-of-concept for the approach, demonstrating it works and scales.
- Built to production quality: tested, documented, reproducible.
- A tool designed to work across chemistries, not just Twist.

## What this project is NOT

- A replacement for all alignment-based pipelines.
- A claim of superiority across all use cases.
- A production-launched commercial tool.

The goal is to show the approach works and is worth using. Not to declare alignment obsolete.
