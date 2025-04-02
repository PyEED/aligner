# Sequence Alignment Tool 🧬

A command-line tool for comparing protein sequences with high performance and flexibility.

## Principle

The tool performs pairwise sequence alignments using either BLOSUM62 or identity scoring.
It supports streaming output and optional pre-filtering based on k-mer matches to improve
performance when dealing with large sequence sets.

## Basic Usage

```bash
./aligner <input> [OPTIONS]
```

## Arguments

| Argument  | Description                                           |
| --------- | ----------------------------------------------------- |
| `<input>` | Path to your input JSON file containing the sequences |

## Options

| Option                    | Description                                                             |
| ------------------------- | ----------------------------------------------------------------------- |
| `-o, --output <FILE>`     | Specify output file path (tab-separated format)                         |
| `-f, --fraction <FLOAT>`  | Set pre-filtering fraction using k-mer matches (0.0-1.0)                |
| `-m, --min-matches <INT>` | Set minimum number of k-mer matches required for alignment (default: 0) |
| `-s, --scoring <TYPE>`    | Choose scoring type: `blosum62` or `identity` (default: identity)       |
| `-h, --help`              | Display help information                                                |
| `-V, --version`           | Show version information                                                |

## Input Format

Your input file should be a JSON file structured as follows:

```json
{
  "Q6A0I3": "MAVMT...",
  "ADV92528.1": "MANPY..."
}
```

## Output Format

The tool generates a tab-separated output with the following columns:

```text
query_id subject_id score seq1_len seq2_len
Q6A0I3 ADV92528.1 ... ... ...
```

## Example Usage

```bash
./aligner input.json -o output.tsv -f 0.5 -s blosum62
```

This will:

- Read sequences from `input.json`
- Use BLOSUM62 scoring matrix
- Apply 50% k-mer pre-filtering
- Save results to `output.tsv`
