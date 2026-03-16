//! Error types for the kam pipeline.
//!
//! All public-facing errors implement [`std::error::Error`] via
//! [`thiserror::Error`]. Library code must never use `unwrap()` — propagate
//! [`KamError`] instead.

use thiserror::Error;

/// Top-level error type for the kam pipeline.
///
/// Individual crates may add their own error variants or wrap `KamError` as
/// needed, but all errors that cross crate boundaries should be expressible as
/// (or convertible to) `KamError`.
///
/// # Example
/// ```
/// use kam_core::error::KamError;
///
/// let err = KamError::ReadTooShort { expected: 10, actual: 3 };
/// assert!(err.to_string().contains("10"));
/// assert!(err.to_string().contains("3"));
/// ```
#[derive(Error, Debug)]
pub enum KamError {
    /// A read was shorter than the minimum required to extract the UMI, skip,
    /// and at least one template base.
    #[error("Read too short: expected at least {expected} bases, got {actual}")]
    ReadTooShort {
        /// Minimum required length.
        expected: usize,
        /// Actual read length.
        actual: usize,
    },

    /// One or more UMI bases fell below the quality threshold and cannot be
    /// used for molecule grouping.
    #[error("Low quality UMI bases")]
    LowQualityUmi,

    /// An underlying I/O error (file not found, permission denied, etc.).
    #[error("IO error: {0}")]
    Io(#[from] std::io::Error),
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn read_too_short_message_contains_lengths() {
        let err = KamError::ReadTooShort { expected: 10, actual: 3 };
        let msg = err.to_string();
        assert!(msg.contains("10"), "expected '10' in: {msg}");
        assert!(msg.contains("3"), "expected '3' in: {msg}");
    }

    #[test]
    fn low_quality_umi_message() {
        let err = KamError::LowQualityUmi;
        assert_eq!(err.to_string(), "Low quality UMI bases");
    }

    #[test]
    fn io_error_wraps_correctly() {
        let io_err = std::io::Error::new(std::io::ErrorKind::NotFound, "file missing");
        let err: KamError = io_err.into();
        assert!(err.to_string().contains("IO error"));
    }
}
