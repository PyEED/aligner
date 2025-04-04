//! Sequence alignment functionality.
//!
//! This module provides functions for performing pairwise sequence alignments,
//! including global alignment and pre-filtering based on k-mer matches.

use bio::alignment::pairwise::*;
use bio::alignment::sparse::find_kmer_matches;
use indicatif::ParallelProgressIterator;
use rayon::ThreadPoolBuilder;
use rayon::prelude::*;
use std::collections::HashMap;
use std::sync::mpsc::Sender;

use crate::utils::setup_progress_bar;

/// Function type for scoring matches between amino acids or nucleotides
pub type MatcherFn = fn(u8, u8) -> i32;

/// Represents the result of a pairwise sequence alignment
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct AlignmentResult {
    /// Identifier of the query sequence
    pub query_id: String,
    /// Identifier of the subject sequence
    pub subject_id: String,
    /// Alignment score, None if alignment was skipped
    pub score: Option<i32>,
    /// Length of sequence 1
    pub seq1_len: usize,
    /// Length of sequence 2
    pub seq2_len: usize,
}

/// Performs pairwise alignments for all unique pairs of sequences in the input,
/// streaming results through a channel.
pub fn align_all_streaming(
    input: &HashMap<String, String>,
    matcher: &MatcherFn,
    fraction: Option<f32>,
    min_matches: usize,
    sender: Sender<AlignmentResult>,
    num_threads: Option<usize>,
) {
    // Set up thread pool if num_threads is specified
    if let Some(n) = num_threads {
        ThreadPoolBuilder::new()
            .num_threads(n)
            .build_global()
            .expect("Failed to initialize thread pool");
    }

    // Use references to keys instead of cloning
    let keys: Vec<&String> = input.keys().collect();

    // Generate all unique pairs using references
    let pairs: Vec<(&String, &String)> = keys
        .iter()
        .enumerate()
        .flat_map(|(i, query_id)| {
            keys[..=i]
                .iter()
                .map(move |subject_id| (*query_id, *subject_id))
        })
        .collect();

    // Setup progress bar with total comparisons
    let progress = setup_progress_bar(pairs.len() as u64);

    // Process alignments in parallel and send results through the channel
    pairs
        .par_iter()
        .progress_with(progress)
        .for_each(|(query_id, subject_id)| {
            if query_id == subject_id {
                return;
            }

            let query_seq = &input[*query_id];
            let subject_seq = &input[*subject_id];
            let score = match fraction {
                Some(fraction) => {
                    if worth_aligning(query_seq, subject_seq, fraction, min_matches) {
                        Some(align(query_seq, subject_seq, matcher))
                    } else {
                        None
                    }
                }
                None => Some(align(query_seq, subject_seq, matcher)),
            };

            let result = AlignmentResult {
                query_id: (*query_id).clone(), // Clone only when creating the result
                subject_id: (*subject_id).clone(), // Clone only when creating the result
                score,
                seq1_len: query_seq.len(),
                seq2_len: subject_seq.len(),
            };

            sender.send(result).expect("Failed to send result");
        });
}

/// Determines if two sequences are worth aligning based on k-mer sharing.
///
/// This function acts as a pre-filter to avoid expensive alignments for sequences
/// that are unlikely to have significant similarity.
///
/// # Arguments
///
/// * `seq1` - First sequence as a string
/// * `seq2` - Second sequence as a string
/// * `fraction` - Fraction of the shorter sequence length to use as k-mer size
///
/// # Returns
///
/// `true` if the sequences share at least one k-mer, `false` otherwise
pub fn worth_aligning(seq1: &str, seq2: &str, fraction: f32, min_matches: usize) -> bool {
    // Use shorter sequence as query
    let (query, subject) = if seq1.len() < seq2.len() {
        (seq1, seq2)
    } else {
        (seq2, seq1)
    };

    let k = (query.len() as f32 * fraction) as usize;
    let kmers = find_kmer_matches(query.as_bytes(), subject.as_bytes(), k);

    kmers.len() >= min_matches
}

/// Performs global alignment between two sequences and returns the alignment score.
///
/// # Arguments
///
/// * `seq1` - First sequence as a string
/// * `seq2` - Second sequence as a string
/// * `matcher` - Scoring function for comparing sequence elements
///
/// # Returns
///
/// The alignment score as an integer
pub fn align(seq1: &str, seq2: &str, matcher: &MatcherFn) -> i32 {
    let mut aligner = Aligner::with_capacity(seq1.len(), seq2.len(), -10, -1, matcher);
    aligner.global(seq1.as_bytes(), seq2.as_bytes()).score
}
