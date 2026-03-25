//! Shared helper for building a [`CallerConfig`] from CLI arguments.
//!
//! Both `call` and `run` subcommands accept the same caller-related flags.
//! This module centralises the construction so neither command silently omits
//! a field.

use kam_call::caller::CallerConfig;

/// Fields shared by `CallArgs` and `RunArgs` that control variant calling.
///
/// Implement this trait for any args struct that carries caller configuration,
/// then call [`caller_config_from_args`] to build the [`CallerConfig`].
pub trait CallerConfigArgs {
    fn min_confidence(&self) -> Option<f64>;
    fn strand_bias_threshold(&self) -> Option<f64>;
    fn min_alt_molecules(&self) -> Option<u32>;
    fn min_alt_duplex(&self) -> Option<u32>;
    fn sv_min_confidence(&self) -> Option<f64>;
    fn sv_min_alt_molecules(&self) -> Option<u32>;
    /// Fisher p-value threshold for strand bias on SV-type variants.
    ///
    /// Unlike the other fields this is not `Option<f64>` — it has a required
    /// default value (1.0, disabled) in the CLI.
    fn sv_strand_bias_threshold(&self) -> f64;
    fn max_vaf(&self) -> Option<f64>;
}

/// Build a [`CallerConfig`] from any args struct that implements
/// [`CallerConfigArgs`].
///
/// Each field falls back to the [`CallerConfig`] default when not supplied by
/// the caller.
pub fn caller_config_from_args(args: &impl CallerConfigArgs) -> CallerConfig {
    let defaults = CallerConfig::default();
    CallerConfig {
        min_confidence: args.min_confidence().unwrap_or(defaults.min_confidence),
        strand_bias_threshold: args
            .strand_bias_threshold()
            .unwrap_or(defaults.strand_bias_threshold),
        min_alt_molecules: args
            .min_alt_molecules()
            .unwrap_or(defaults.min_alt_molecules),
        min_alt_duplex: args.min_alt_duplex().unwrap_or(defaults.min_alt_duplex),
        sv_min_confidence: args
            .sv_min_confidence()
            .unwrap_or(defaults.sv_min_confidence),
        sv_min_alt_molecules: args
            .sv_min_alt_molecules()
            .unwrap_or(defaults.sv_min_alt_molecules),
        sv_strand_bias_threshold: args.sv_strand_bias_threshold(),
        max_vaf: args.max_vaf().or(defaults.max_vaf),
        ..defaults
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    struct TestArgs {
        min_confidence: Option<f64>,
        strand_bias_threshold: Option<f64>,
        min_alt_molecules: Option<u32>,
        min_alt_duplex: Option<u32>,
        sv_min_confidence: Option<f64>,
        sv_min_alt_molecules: Option<u32>,
        sv_strand_bias_threshold: f64,
        max_vaf: Option<f64>,
    }

    impl CallerConfigArgs for TestArgs {
        fn min_confidence(&self) -> Option<f64> {
            self.min_confidence
        }
        fn strand_bias_threshold(&self) -> Option<f64> {
            self.strand_bias_threshold
        }
        fn min_alt_molecules(&self) -> Option<u32> {
            self.min_alt_molecules
        }
        fn min_alt_duplex(&self) -> Option<u32> {
            self.min_alt_duplex
        }
        fn sv_min_confidence(&self) -> Option<f64> {
            self.sv_min_confidence
        }
        fn sv_min_alt_molecules(&self) -> Option<u32> {
            self.sv_min_alt_molecules
        }
        fn sv_strand_bias_threshold(&self) -> f64 {
            self.sv_strand_bias_threshold
        }
        fn max_vaf(&self) -> Option<f64> {
            self.max_vaf
        }
    }

    fn all_none_args() -> TestArgs {
        TestArgs {
            min_confidence: None,
            strand_bias_threshold: None,
            min_alt_molecules: None,
            min_alt_duplex: None,
            sv_min_confidence: None,
            sv_min_alt_molecules: None,
            sv_strand_bias_threshold: 1.0,
            max_vaf: None,
        }
    }

    /// All-None args produce a config equal to CallerConfig::default().
    #[test]
    fn all_none_produces_defaults() {
        let args = all_none_args();
        let cfg = caller_config_from_args(&args);
        let defaults = CallerConfig::default();
        assert_eq!(cfg.min_confidence, defaults.min_confidence);
        assert_eq!(cfg.strand_bias_threshold, defaults.strand_bias_threshold);
        assert_eq!(cfg.min_alt_molecules, defaults.min_alt_molecules);
        assert_eq!(cfg.min_alt_duplex, defaults.min_alt_duplex);
        assert_eq!(cfg.sv_min_confidence, defaults.sv_min_confidence);
        assert_eq!(cfg.sv_min_alt_molecules, defaults.sv_min_alt_molecules);
        assert_eq!(cfg.max_vaf, defaults.max_vaf);
    }

    /// Provided values override defaults.
    #[test]
    fn provided_values_override_defaults() {
        let args = TestArgs {
            min_confidence: Some(0.99),
            strand_bias_threshold: Some(0.001),
            min_alt_molecules: Some(5),
            min_alt_duplex: Some(2),
            sv_min_confidence: Some(0.80),
            sv_min_alt_molecules: Some(3),
            sv_strand_bias_threshold: 0.05,
            max_vaf: Some(0.35),
        };
        let cfg = caller_config_from_args(&args);
        assert_eq!(cfg.min_confidence, 0.99);
        assert_eq!(cfg.strand_bias_threshold, 0.001);
        assert_eq!(cfg.min_alt_molecules, 5);
        assert_eq!(cfg.min_alt_duplex, 2);
        assert_eq!(cfg.sv_min_confidence, 0.80);
        assert_eq!(cfg.sv_min_alt_molecules, 3);
        assert_eq!(cfg.sv_strand_bias_threshold, 0.05);
        assert_eq!(cfg.max_vaf, Some(0.35));
    }
}
