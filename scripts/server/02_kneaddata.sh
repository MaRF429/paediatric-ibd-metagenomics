#!/usr/bin/env bash
set -euo pipefail

# KneadData v0.12.2; process one sample with an hg38 Bowtie2 database.
# Usage: bash 02_kneaddata.sh paired READ1 READ2 OUTPUT_DIR HG38_DB [THREADS] [PROCESSES]
# Usage: bash 02_kneaddata.sh single FASTQ OUTPUT_DIR HG38_DB [THREADS] [PROCESSES]

run_paired() {
read1="${1:?Read 1 FASTQ is required}"
read2="${2:?Read 2 FASTQ is required}"
output_dir="${3:?Output directory is required}"
hg38_db="${4:?hg38 Bowtie2 database path is required}"
threads="${5:-24}"
processes="${6:-8}"
[[ -r "$read1" && -r "$read2" ]] || { echo "FASTQ inputs must be readable." >&2; exit 1; }
[[ "$threads" =~ ^[1-9][0-9]*$ && "$processes" =~ ^[1-9][0-9]*$ ]] || { echo "Threads and processes must be positive." >&2; exit 1; }

# TruSeq.fa must be available to Trimmomatic as in the original workflow.
mkdir -p "$output_dir"
kneaddata \
    --input1 "$read1" \
    --input2 "$read2" \
    --output "$output_dir" \
    --reference-db "$hg38_db" \
    --verbose \
    --remove-intermediate-output \
    --threads "$threads" \
    --processes "$processes" \
    --trimmomatic-options "ILLUMINACLIP:TruSeq.fa:2:30:10:2:TRUE LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 HEADCROP:8 MINLEN:36" \
    > "$output_dir/kneaddata.log" 2>&1
}

run_single() {
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
}

mode="${1:?Mode must be paired or single}"
shift
case "$mode" in
    paired) run_paired "$@" ;;
    single) run_single "$@" ;;
    *) echo "Mode must be paired or single." >&2; exit 1 ;;
esac
