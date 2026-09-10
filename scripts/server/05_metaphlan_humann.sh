#!/usr/bin/env bash
set -euo pipefail

# MetaPhlAn v4.0.6 and HUMAnN v3.9; one sample per invocation.
# Usage: bash 05_metaphlan_humann.sh SAMPLE_ID COHORT FASTQ OUTPUT_ROOT METAPHLAN_DB CHOCOPHLAN_DB UNIREF90_DB [THREADS]
sample_id="${1:?Sample ID is required}"
cohort="${2:?Cohort accession is required}"
input_fastq="${3:?Merged FASTQ is required}"
output_root="${4:?Output root is required}"
metaphlan_db="${5:?MetaPhlAn database directory is required}"
chocophlan_db="${6:?ChocoPhlAn database directory is required}"
uniref90_db="${7:?UniRef90 database directory is required}"
threads="${8:-20}"
[[ "$sample_id" =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ ]] || { echo "Invalid sample ID." >&2; exit 1; }
[[ -r "$input_fastq" ]] || { echo "FASTQ input must be readable." >&2; exit 1; }
[[ -d "$metaphlan_db" && -d "$chocophlan_db" && -d "$uniref90_db" ]] || { echo "Database directories must exist." >&2; exit 1; }
[[ "$threads" =~ ^[1-9][0-9]*$ ]] || { echo "Threads must be positive." >&2; exit 1; }

case "$cohort" in
    PRJNA398089|PRJNA1265906|PRJNA1045596|PRJNA389280|HRA007915)
        group="discovery" ;;
    PRJEB76677|SRP057027|PRJNA759642|PRJNA922068)
        group="independent_validation" ;;
    *) echo "Unknown cohort accession; no group assigned." >&2; exit 1 ;;
esac

output_dir="$output_root/$group/$cohort/$sample_id"
bugs="$output_dir/${sample_id}_metaphlan_bugs.tsv"
mkdir -p "$output_dir"
metaphlan "$input_fastq" \
    --input_type fastq \
    --bowtie2db "$metaphlan_db" \
    --index mpa_vJun23_CHOCOPhlAnSGB_202403 \
    --nproc "$threads" \
    --offline \
    --bowtie2out "$output_dir/${sample_id}.bowtie2out.txt" \
    -o "$bugs" \
    2> "$output_dir/${sample_id}_metaphlan.stderr.log"
[[ -s "$bugs" ]] || { echo "MetaPhlAn profile is empty." >&2; exit 1; }

humann \
    --input "$input_fastq" \
    --taxonomic-profile "$bugs" \
    --output "$output_dir" \
    --output-basename "$sample_id" \
    --threads "$threads" \
    --remove-temp-output \
    --nucleotide-database "$chocophlan_db" \
    --protein-database "$uniref90_db" \
    --o-log "$output_dir/${sample_id}.humann.log" \
    --log-level INFO
