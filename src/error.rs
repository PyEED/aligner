//! Error handling module for the sequence aligner application.
//!
//! This module defines the `AlignerError` type, which represents all possible
//! errors that can occur during the execution of the application. It implements
//! the `Error` trait through the `thiserror` crate's derive macro.

/// Error type for the aligner application.
///
/// This enum represents all possible errors that can occur in the application.
/// It implements the standard Error trait through thiserror's derive macro.
#[derive(Debug, thiserror::Error)]
pub(crate) enum AlignerError {
    /// IO error that occurs during file operations.
    ///
    /// This variant wraps a standard IO error and is typically returned
    /// when there are issues opening or reading input files.
    #[error("IO error: {0}")]
    Io(#[from] std::io::Error),

    /// Parse error that occurs during JSON deserialization.
    ///
    /// This variant wraps a serde_json error and is returned when
    /// the input JSON file cannot be properly parsed into the expected format.
    #[error("Parse error: {0}")]
    Parse(#[from] serde_json::Error),
}
