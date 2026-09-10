#!/usr/bin/env bash
set -euo pipefail

# Analyze one input table per invocation.
# Usage: bash 07_fastspar.sh INPUT_TSV OUTPUT_DIR PREFIX [THREADS]
input_table="${1:?Input TSV table is required}"
output_dir="${2:?Output directory is required}"
prefix="${3:?Output prefix is required}"
threads="${4:-4}"
[[ -f "$input_table" && -r "$input_table" ]] || { echo "Input table must be readable." >&2; exit 1; }
[[ "$prefix" =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ ]] || { echo "Invalid output prefix." >&2; exit 1; }
[[ "$threads" =~ ^[1-9][0-9]*$ ]] || { echo "Threads must be positive." >&2; exit 1; }

mkdir -p "$output_dir/bootstrap_counts" "$output_dir/bootstrap_correlation"
fastspar \
    --otu_table "$input_table" \
    --iterations 50 \
    --threads "$threads" \
    --correlation "$output_dir/${prefix}_correlation.tsv" \
    --covariance "$output_dir/${prefix}_covariance.tsv"

fastspar_bootstrap \
    --otu_table "$input_table" \
    --number 1000 \
    --prefix "$output_dir/bootstrap_counts/$prefix"

# Quote command arguments when GNU Parallel invokes its child shell.
shopt -s failglob
parallel --quote fastspar \
    --otu_table '{}' \
    --iterations 5 \
    --threads "$threads" \
    --correlation "$output_dir/bootstrap_correlation/cor_{/}" \
    --covariance "$output_dir/bootstrap_correlation/cov_{/}" \
    ::: "$output_dir/bootstrap_counts/${prefix}_"*.tsv

# Matches cor_<PREFIX>_<N>.tsv generated above.
fastspar_pvalues \
    --otu_table "$input_table" \
    --correlation "$output_dir/${prefix}_correlation.tsv" \
    --prefix "$output_dir/bootstrap_correlation/cor_${prefix}_" \
    --permutations 1000 \
    --outfile "$output_dir/${prefix}_pvalues.tsv"
