#!/usr/bin/env bash
set -euo pipefail

# MetaPhlAn v4.0.6 and HUMAnN v3.9.
# Profile one sample or build tables separately within one group.
# Usage: bash 04_metaphlan_humann.sh profile SAMPLE_ID COHORT FASTQ OUTPUT_ROOT METAPHLAN_DB CHOCOPHLAN_DB UNIREF90_DB [THREADS]
# Usage: bash 04_metaphlan_humann.sh tables GROUP PROFILE_ROOT TABLE_ROOT
# GROUP is discovery or independent_validation.
# Profile mode expects an already prepared FASTQ; paired reads are not merged here.
# Independent validation never informs discovery table construction.

run_profile() {
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
}

build_tables() {
group="${1:?Group is required}"
profile_root="${2:?Profile root is required}"
table_root="${3:?Table output root is required}"
case "$group" in
    discovery|independent_validation) ;;
    *) echo "Group must be discovery or independent_validation." >&2; exit 1 ;;
esac
input_dir="$profile_root/$group"
output_dir="$table_root/$group"
[[ -d "$input_dir" ]] || { echo "Group profile directory does not exist." >&2; exit 1; }
shopt -s failglob
mkdir -p "$output_dir"

humann_join_tables --input "$input_dir" --file_name genefamilies \
    --output "$output_dir/${group}_genefamilies.tsv" --search-subdirectories
humann_join_tables --input "$input_dir" --file_name pathabundance \
    --output "$output_dir/${group}_pathabundance.tsv" --search-subdirectories

humann_regroup_table --input "$output_dir/${group}_genefamilies.tsv" \
    --groups uniref90_ko --output "$output_dir/${group}_KO.tsv"
humann_regroup_table --input "$output_dir/${group}_genefamilies.tsv" \
    --groups uniref90_eggnog --output "$output_dir/${group}_EggNOG.tsv"

humann_renorm_table --input "$output_dir/${group}_pathabundance.tsv" \
    --units relab --output "$output_dir/${group}_pathway_relab.tsv"
humann_renorm_table --input "$output_dir/${group}_KO.tsv" \
    --units relab --output "$output_dir/${group}_KO_relab.tsv"
humann_renorm_table --input "$output_dir/${group}_EggNOG.tsv" \
    --units relab --output "$output_dir/${group}_EggNOG_relab.tsv"

humann_split_stratified_table --input "$output_dir/${group}_pathway_relab.tsv" \
    --output "$output_dir/${group}_pathway_tables"
humann_split_stratified_table --input "$output_dir/${group}_KO_relab.tsv" \
    --output "$output_dir/${group}_KO_tables"
humann_split_stratified_table --input "$output_dir/${group}_EggNOG_relab.tsv" \
    --output "$output_dir/${group}_EggNOG_tables"

awk 'NR==1 || $1 !~ /READS_UNMAPPED|UNMAPPED|UNGROUPED|unknown/' \
    "$output_dir/${group}_KO_tables/${group}_KO_relab_unstratified.tsv" \
    > "$output_dir/${group}_KO_clean.tsv"
awk 'NR==1 || $1 !~ /READS_UNMAPPED|UNMAPPED|UNGROUPED|unknown/' \
    "$output_dir/${group}_EggNOG_tables/${group}_EggNOG_relab_unstratified.tsv" \
    > "$output_dir/${group}_EggNOG_clean.tsv"
awk 'NR==1 || $1 !~ /READS_UNMAPPED|UNMAPPED|UNINTEGRATED|UNGROUPED/' \
    "$output_dir/${group}_pathway_tables/${group}_pathway_relab_unstratified.tsv" \
    > "$output_dir/${group}_pathway_clean.tsv"

# Expand exactly the cohort/sample directory levels within this group.
merge_metaphlan_tables.py "$input_dir"/*/*/*_metaphlan_bugs.tsv \
    > "$output_dir/${group}_metaphlan.tsv"
}

mode="${1:?Mode must be profile or tables}"
shift
case "$mode" in
    profile) run_profile "$@" ;;
    tables) build_tables "$@" ;;
    *) echo "Mode must be profile or tables." >&2; exit 1 ;;
esac
