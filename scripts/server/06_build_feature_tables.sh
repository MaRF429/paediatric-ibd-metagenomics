#!/usr/bin/env bash
set -euo pipefail

# Build tables separately; independent validation never informs discovery.
# Usage: bash 06_build_feature_tables.sh GROUP PROFILE_ROOT TABLE_ROOT
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
