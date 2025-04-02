# Basic usage
./aligner test_input.json

# With output file
./aligner test_input.json -o output.tsv

# With BLOSUM62 and fraction 0.035
./aligner test_input.json -o output_blosum62.tsv -f 0.035 -s blosum62

# With identity scoring and fraction 0.035
./aligner test_input.json -o output_identity.tsv -f 0.035 -s identity
