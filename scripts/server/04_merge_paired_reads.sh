#!/usr/bin/env bash
set -euo pipefail

# Concatenate uncompressed KneadData paired FASTQ outputs; retain both inputs.
# Usage: bash 04_merge_paired_reads.sh READ1_FASTQ READ2_FASTQ SAMPLE_ID OUTPUT_DIR
read1="${1:?Clean paired Read 1 is required}"
read2="${2:?Clean paired Read 2 is required}"
sample_id="${3:?Sample ID is required}"
output_dir="${4:?Output directory is required}"
[[ -r "$read1" && -r "$read2" ]] || { echo "FASTQ inputs must be readable." >&2; exit 1; }
case "$read1:$read2" in
    *.gz:*|*:*.gz) echo "Inputs must be uncompressed FASTQ files." >&2; exit 1 ;;
esac
[[ "$sample_id" =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ ]] || { echo "Invalid sample ID." >&2; exit 1; }

mkdir -p "$output_dir"
cat -- "$read1" "$read2" | gzip -c > "$output_dir/${sample_id}_kneaddata_paired.fastq.gz"
