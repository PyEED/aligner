//! A command-line sequence alignment tool for comparing protein or nucleotide sequences.
//!
//! This tool performs pairwise sequence alignments using either BLOSUM62 or identity scoring.
//! It supports streaming output and optional pre-filtering based on k-mer matches to improve
//! performance when dealing with large sequence sets.
//!
//! # Usage
//!
//! ```text
//! aligner <input> [OPTIONS]
//!
//! Arguments:
//!   <input>    Path to input JSON file containing sequences
//!
//! Options:
//!   -o, --output <FILE>     Path to output file (tab-separated format)
//!   -f, --fraction <FLOAT>  Fraction for pre-filtering using k-mer matches (0.0-1.0)
//!   -s, --scoring <TYPE>    Scoring type: blosum62 or identity [default: identity]
//!   -h, --help             Print help
//!   -V, --version          Print version
//! ```
//!
//! # Example
//!
//! ```bash
//! aligner input.json -o output.tsv -f 0.5 -s blosum62
//! ```
//!
//! # Input Format
//!
//! The input JSON file should have the following format:
//!
//! ```json
//! {
//!   "Q6A0I3": "MAVMT...",
//!   "ADV92528.1": "MANPY..."
//! }
//! ```
//!
//! # Output Format
//!
//! The output file will be tab-separated with the following columns:
//!
//! ```text
//! query_id\tsubject_id\tscore\tseq1_len\tseq2_len
//! Q6A0I3\tADV92528.1\t...\t...\t...
//! ```

mod align;
mod error;
mod utils;

use align::{MatcherFn, align_all_streaming};
use bio::scores::blosum62;
use clap::{Parser, ValueEnum};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;
use std::sync::mpsc;
use std::time::Instant;
use utils::parse_input;

/// Supported scoring matrices for sequence alignment
#[derive(Debug, Copy, Clone, ValueEnum)]
enum ScoringType {
    /// BLOSUM62 scoring matrix, optimized for protein sequence alignments
    Blosum62,
    /// Simple identity scoring: match=1, mismatch=-1
    Identity,
}

/// Command-line arguments for the sequence alignment tool
#[derive(Parser, Debug)]
#[command(author, version, about = "Sequence alignment tool")]
struct Args {
    /// Path to input JSON file containing sequences.
    /// The file should contain a JSON object where keys are sequence identifiers
    /// and values are the sequences as strings.
    #[arg(help = "Path to input JSON file containing sequences")]
    input: PathBuf,

    /// Path to output file (optional).
    /// If provided, results will be written in tab-separated format with columns:
    /// query_id, subject_id, score, seq1_len, seq2_len
    #[arg(short, long, help = "Path to output file")]
    output: Option<PathBuf>,

    /// Fraction for pre-filtering sequences using k-mer matches (between 0 and 1).
    /// Higher values are more stringent. If provided, sequences sharing fewer k-mers
    /// than this threshold will be skipped, improving performance.
    #[arg(short, long, help = "Fraction for pre-filtering using k-mer matches")]
    fraction: Option<f32>,

    /// Scoring type to use for alignment.
    /// BLOSUM62 is recommended for protein sequences, while Identity scoring
    /// works for both protein and nucleotide sequences.
    #[arg(short, long, value_enum, default_value_t = ScoringType::Identity, help = "Scoring type to use for alignment")]
    scoring: ScoringType,

    /// Minimum number of k-mer matches required for alignment.
    #[arg(
        short,
        long,
        default_value = "0",
        help = "Minimum number of k-mer matches required for alignment"
    )]
    min_matches: usize,

    /// Number of threads to use for parallel processing.
    #[arg(
        short,
        long,
        help = "Number of threads to use for parallel processing. If not provided, the number of threads will be determined automatically."
    )]
    threads: Option<usize>,
}

/// Scoring function wrapper that supports built-in and custom scoring matrices
pub enum Matcher {
    /// BLOSUM62 scoring matrix from the bio crate
    Blosum62,
    /// Simple identity scoring (match=1, mismatch=-1)
    Identity,
    /// Custom scoring function
    Custom(MatcherFn),
}

impl Matcher {
    /// Returns the actual scoring function for use in alignment
    fn score(&self) -> MatcherFn {
        match self {
            Matcher::Blosum62 => blosum62,
            Matcher::Identity => |a, b| if a == b { 1 } else { -1 },
            Matcher::Custom(matcher) => *matcher,
        }
    }
}

fn main() {
    let args = Args::parse();

    // Validate fraction if provided
    if let Some(fraction) = args.fraction {
        if !(0.0..=1.0).contains(&fraction) {
            eprintln!("Error: fraction must be between 0 and 1");
            std::process::exit(1);
        }
    }

    let input = match parse_input(&args.input) {
        Ok(input) => input,
        Err(e) => {
            eprintln!("Error reading input file: {}", e);
            std::process::exit(1);
        }
    };

    let match_fn = match args.scoring {
        ScoringType::Blosum62 => Matcher::Blosum62.score(),
        ScoringType::Identity => Matcher::Identity.score(),
    };

    let start = Instant::now();

    // Set up output writer if path is specified
    let mut writer = args.output.map(|path| {
        let file = File::create(path).expect("Failed to create output file");
        let mut writer = BufWriter::new(file);
        // Write CSV header
        writeln!(writer, "query_id\tsubject_id\tscore\tseq1_len\tseq2_len")
            .expect("Failed to write header");
        writer
    });

    // Create channel for streaming results
    let (tx, rx) = mpsc::channel();

    // Spawn the alignment computation using rayon's threading
    let computation_handle = std::thread::spawn(move || {
        align_all_streaming(
            &input,
            &match_fn,
            args.fraction,
            args.min_matches,
            tx,
            args.threads,
        )
    });

    // Process results as they arrive
    let mut total_results = 0;
    for result in rx {
        total_results += 1;
        if let Some(ref mut w) = writer {
            writeln!(
                w,
                "{}\t{}\t{}\t{}\t{}",
                result.query_id,
                result.subject_id,
                result.score.unwrap_or(-1),
                result.seq1_len,
                result.seq2_len
            )
            .expect("Failed to write result");
        }
    }

    // Wait for computation to finish
    computation_handle
        .join()
        .expect("Computation thread panicked");

    let duration = start.elapsed().as_secs_f32();
    println!("Processed {} alignments in {:.2}s", total_results, duration);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_input() {
        let input = parse_input("tests/data/test_input.json").unwrap();
        assert_eq!(input.len(), 2);
        assert!(input.contains_key("Q6A0I3"));
        assert!(input.contains_key("ADV92528.1"));
        assert_eq!(
            input["Q6A0I3"],
            "MAVMTPRRERSSLLSRALRFTAAAATALVTAVSLAAPAHAANPYERGPNPTDALLEARSGPFSVSEERASRFGADGFGGGTIYYPRENNTYGAVAISPGYTGTQASVAWLGKRIASHGFVVITIDTNTTLDQPDSRARQLNAALDYMINDASSAVRSRIDSSRLAVMGHSMGGGGSLRLASQRPDLKAAIPLTPWHLNKNWSSVRVPTLIIGADLDTIAPVLTHARPFYNSLPTSISKAYLELDGATHFAPNIPNKIIGKYSVAWLKRFVDNDTRYTQFLCPGPRDGLFGEVEEYRSTCPF"
        );
        assert_eq!(
            input["ADV92528.1"],
            "MANPYERGPNPTDALLEARSGPFSVSEENVSRLSASGFGGGTIYYPRENNTYGAVAISPGYTGTEASIAWLGERIASHGFVVITIDTITTLDQPDSRAEQLNAALNHMINRASSTVRSRIDSSRLAVMGHSMGGGGSLRLASQRPDLKAAIPLTPWHLNKNWSSVRVPTLIIGADLDTIAPVLTHARPFYNSLPTSISKAYLELDGATHFAPNIPNKIIGKYSVAWLKRFVDNDTRYTQFLCPGPRDGLFGEVEEYRSTCPF"
        );
    }
}
