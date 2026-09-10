#!/usr/bin/env bash
set -euo pipefail

# KneadData v0.12.2; one single-end sample with an hg38 Bowtie2 database.
# Usage: bash 03_kneaddata_single.sh FASTQ OUTPUT_DIR HG38_DB [THREADS] [PROCESSES]
read="${1:?FASTQ input is required}"
output_dir="${2:?Output directory is required}"
hg38_db="${3:?hg38 Bowtie2 database path is required}"
threads="${4:-24}"
processes="${5:-8}"
[[ -r "$read" ]] || { echo "FASTQ input must be readable." >&2; exit 1; }
[[ "$threads" =~ ^[1-9][0-9]*$ && "$processes" =~ ^[1-9][0-9]*$ ]] || { echo "Threads and processes must be positive." >&2; exit 1; }

# TruSeq.fa must be available to Trimmomatic as in the original workflow.
mkdir -p "$output_dir"
kneaddata \
    --unpaired "$read" \
    --output "$output_dir" \
    --reference-db "$hg38_db" \
    --verbose \
    --remove-intermediate-output \
    --threads "$threads" \
    --processes "$processes" \
    --trimmomatic-options "ILLUMINACLIP:TruSeq.fa:2:30:10:2:TRUE LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 HEADCROP:8 MINLEN:36" \
    > "$output_dir/kneaddata.log" 2>&1
