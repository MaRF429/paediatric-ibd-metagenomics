#!/usr/bin/env bash
set -euo pipefail

# FastQC v0.12.1 and MultiQC v1.29.
# Usage: bash 01_fastqc_multiqc.sh FASTQ_DIR OUTPUT_DIR [THREADS]
input_dir="${1:?FASTQ input directory is required}"
output_dir="${2:?Output directory is required}"
threads="${3:-4}"
[[ -d "$input_dir" ]] || { echo "FASTQ directory does not exist." >&2; exit 1; }
[[ "$threads" =~ ^[1-9][0-9]*$ ]] || { echo "Threads must be positive." >&2; exit 1; }

# Reject an unmatched input glob instead of passing it to FastQC.
shopt -s failglob
mkdir -p "$output_dir"
fastqc "$input_dir"/*.fastq.gz -t "$threads" -o "$output_dir"
multiqc "$output_dir" -o "$output_dir"
