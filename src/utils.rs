//! Utility functions for the sequence aligner application.
//!
//! This module provides helper functions for progress tracking and input parsing.

use indicatif::ProgressBar;
use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;
use std::path::PathBuf;

use crate::error::AlignerError;

/// Creates and configures a progress bar for tracking alignment operations.
///
/// This function sets up a progress bar with a custom style to display the
/// progress of sequence alignments in a user-friendly format.
///
/// # Arguments
///
/// * `total_comparisons` - The total number of comparisons to be performed
///
/// # Returns
///
/// A configured `ProgressBar` instance ready for tracking progress
pub fn setup_progress_bar(total_comparisons: u64) -> ProgressBar {
    let progress = ProgressBar::new(total_comparisons);
    progress.set_style(
        indicatif::ProgressStyle::default_bar()
            .template(
                "{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} ({eta})",
            )
            .unwrap()
            .progress_chars("#>-"),
    );
    progress
}

/// Parses a JSON input file containing sequence data.
///
/// Reads a JSON file where keys are sequence identifiers and values are the
/// actual sequences, and converts it into a HashMap for efficient lookup.
///
/// # Arguments
///
/// * `path` - Path to the JSON file containing sequence data
///
/// # Returns
///
/// * `Result<HashMap<String, String>, AlignerError>` - A map of sequence IDs to sequences,
///   or an error if the file cannot be opened or parsed
///
/// # Errors
///
/// Returns `AlignerError::Io` if the file cannot be opened or read.
/// Returns `AlignerError::Parse` if the JSON is malformed or doesn't match the expected format.
pub fn parse_input(path: impl Into<PathBuf>) -> Result<HashMap<String, String>, AlignerError> {
    let content = File::open(path.into()).map_err(AlignerError::Io)?;
    let reader = BufReader::new(content);
    serde_json::from_reader(reader).map_err(AlignerError::Parse)
}
