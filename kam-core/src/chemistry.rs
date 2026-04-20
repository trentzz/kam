//! Chemistry-specific read structure definitions.
//!
//! `kam` currently targets Twist Bioscience duplex UMI chemistry (`5M2S+T`),
//! but [`ReadStructure`] is kept general enough to accommodate other
//! UMI-based chemistries in the future.

/// Describes how a read is structured: UMI bases, skip/spacer bases, and the
/// remaining template sequence.
///
/// For Twist duplex UMI chemistry the read layout on both R1 and R2 is:
/// ```text
/// [5 bp UMI][2 bp skip][...template...]
/// ```
///
/// Use [`ReadStructure::twist_umi_duplex`] to obtain the preset for this
/// chemistry, or construct a custom `ReadStructure` directly.
///
/// # Example
/// ```
/// use kam_core::chemistry::ReadStructure;
///
/// let rs = ReadStructure::twist_umi_duplex();
/// assert_eq!(rs.template_start(), 7);
/// ```
#[derive(Debug, Clone, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub struct ReadStructure {
    /// Length of the UMI in bases.
    pub umi_length: usize,
    /// Length of the skip (spacer) region immediately after the UMI.
    /// The template begins after the skip.
    pub skip_length: usize,
}

impl ReadStructure {
    /// Return the [`ReadStructure`] preset for Twist Bioscience duplex UMI
    /// chemistry: 5 bp UMI followed by a 2 bp skip (`5M2S+T`).
    ///
    /// # Example
    /// ```
    /// use kam_core::chemistry::ReadStructure;
    ///
    /// let rs = ReadStructure::twist_umi_duplex();
    /// assert_eq!(rs.umi_length, 5);
    /// assert_eq!(rs.skip_length, 2);
    /// ```
    pub fn twist_umi_duplex() -> Self {
        Self {
            umi_length: 5,
            skip_length: 2,
        }
    }

    /// Return the [`ReadStructure`] preset for simplex 12 bp UMI chemistry:
    /// 12 bp UMI with no skip region (`12M+T`).
    ///
    /// # Example
    /// ```
    /// use kam_core::chemistry::ReadStructure;
    ///
    /// let rs = ReadStructure::simplex_12bp();
    /// assert_eq!(rs.umi_length, 12);
    /// assert_eq!(rs.skip_length, 0);
    /// assert_eq!(rs.template_start(), 12);
    /// ```
    pub fn simplex_12bp() -> Self {
        Self {
            umi_length: 12,
            skip_length: 0,
        }
    }

    /// Return the byte offset at which the template sequence begins within a
    /// raw read.
    ///
    /// This is simply `umi_length + skip_length`.
    ///
    /// # Example
    /// ```
    /// use kam_core::chemistry::ReadStructure;
    ///
    /// assert_eq!(ReadStructure::twist_umi_duplex().template_start(), 7);
    /// ```
    pub fn template_start(&self) -> usize {
        self.umi_length + self.skip_length
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn twist_umi_duplex_umi_length() {
        assert_eq!(ReadStructure::twist_umi_duplex().umi_length, 5);
    }

    #[test]
    fn twist_umi_duplex_skip_length() {
        assert_eq!(ReadStructure::twist_umi_duplex().skip_length, 2);
    }

    #[test]
    fn template_start_is_7_for_twist() {
        assert_eq!(ReadStructure::twist_umi_duplex().template_start(), 7);
    }

    #[test]
    fn custom_read_structure_template_start() {
        let rs = ReadStructure {
            umi_length: 8,
            skip_length: 3,
        };
        assert_eq!(rs.template_start(), 11);
    }

    #[test]
    fn simplex_12bp_umi_length() {
        assert_eq!(ReadStructure::simplex_12bp().umi_length, 12);
    }

    #[test]
    fn simplex_12bp_skip_length() {
        assert_eq!(ReadStructure::simplex_12bp().skip_length, 0);
    }

    #[test]
    fn simplex_12bp_template_start_is_12() {
        assert_eq!(ReadStructure::simplex_12bp().template_start(), 12);
    }
}
