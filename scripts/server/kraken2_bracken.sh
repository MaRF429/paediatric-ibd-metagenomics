#!/usr/bin/env bash
set -euo pipefail

# Classify one sample, estimate species abundance, and convert its report.
# Usage: bash kraken2_bracken.sh MODE READ1 READ2 SAMPLE_ID DB_DIR OUTPUT_DIR [THREADS]
# MODE is paired or single; use - as READ2 in single mode.
# Supply an existing database with the matching Bracken distribution files.
mode="${1:?Mode is required}"
read1="${2:?Read 1 FASTQ is required}"
read2="${3:?Read 2 FASTQ or - is required}"
sample_id="${4:?Sample ID is required}"
db_dir="${5:?Database directory is required}"
output_dir="${6:?Output directory is required}"
threads="${7:-4}"
[[ -f "$read1" && -r "$read1" ]] || { echo "Read 1 must be readable." >&2; exit 1; }
[[ -d "$db_dir" ]] || { echo "Database directory does not exist." >&2; exit 1; }
[[ "$sample_id" =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ ]] || { echo "Invalid sample ID." >&2; exit 1; }
[[ "$threads" =~ ^[1-9][0-9]*$ ]] || { echo "Threads must be positive." >&2; exit 1; }
case "$mode" in
    paired)
        [[ -f "$read2" && -r "$read2" ]] || { echo "Read 2 must be readable." >&2; exit 1; } ;;
    single)
        [[ "$read2" == "-" ]] || { echo "Use - as Read 2 in single mode." >&2; exit 1; } ;;
    *) echo "Mode must be paired or single." >&2; exit 1 ;;
esac

mkdir -p "$output_dir/kraken2" "$output_dir/bracken"
case "$mode" in
    paired)
        kraken2 --db "$db_dir" --paired --threads "$threads" \
            --report "$output_dir/kraken2/${sample_id}.kreport" \
            --output "$output_dir/kraken2/${sample_id}.kraken" \
            --use-names "$read1" "$read2" ;;
    single)
        kraken2 --db "$db_dir" --threads "$threads" \
            --report "$output_dir/kraken2/${sample_id}.kreport" \
            --output "$output_dir/kraken2/${sample_id}.kraken" \
            --use-names "$read1" ;;
esac

# Preserve the species level and threshold; -t is not a thread count.
# The original command leaves the Bracken read-length setting at its default.
bracken -d "$db_dir" \
    -i "$output_dir/kraken2/${sample_id}.kreport" \
    -o "$output_dir/bracken/${sample_id}.bracken.S" \
    -w "$output_dir/bracken/${sample_id}.bracken.S.kreport" \
    -l S -t 12

# Convert only this sample's report; do not combine unrelated samples.
kraken-biom "$output_dir/bracken/${sample_id}.bracken.S.kreport" \
    --max D --output_fp "$output_dir/bracken/${sample_id}.S.biom"
biom convert -i "$output_dir/bracken/${sample_id}.S.biom" \
    -o "$output_dir/bracken/${sample_id}.S.tsv" --to-tsv
biom convert -i "$output_dir/bracken/${sample_id}.S.biom" \
    -o "$output_dir/bracken/${sample_id}.S.json" --to-json
