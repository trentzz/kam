//! Statistical variant calling with confidence intervals.
//!
//! `kam-call` applies binomial and Beta-posterior models to molecule-level
//! evidence from the de Bruijn graph to produce calibrated variant calls
//! with posterior credible intervals and strand-bias filters.

pub mod allele;
pub mod caller;
pub mod fusion;
pub mod output;
pub mod targeting;
