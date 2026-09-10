#!/usr/bin/env bash
set -euo pipefail

# Run the six specified species-function association analyses.
# Usage: bash 08_halla.sh INPUT_DIR OUTPUT_ROOT
input_dir="${1:?Input table directory is required}"
output_root="${2:?Output root is required}"
[[ -d "$input_dir" ]] || { echo "Input directory does not exist." >&2; exit 1; }
mkdir -p "$output_root"

halla -x "$input_dir/IBD_species.tsv" -y "$input_dir/IBD_KO.tsv" \
    -o "$output_root/halla_KO_IBD" \
    -m spearman --fdr_alpha 0.05 --permute_iters 1000

halla -x "$input_dir/IBD_species.tsv" -y "$input_dir/IBD_EggNOG.tsv" \
    -o "$output_root/halla_EggNOG_IBD" \
    -m spearman --fdr_alpha 0.05 --permute_iters 1000

halla -x "$input_dir/IBD_species.tsv" -y "$input_dir/IBD_pathway.tsv" \
    -o "$output_root/halla_pathway_IBD" \
    -m spearman --fdr_alpha 0.05 --permute_iters 1000

halla -x "$input_dir/CDUC_species.tsv" -y "$input_dir/CDUC_KO.tsv" \
    -o "$output_root/halla_KO_CDUC" \
    -m spearman --fdr_alpha 0.05 --permute_iters 1000

halla -x "$input_dir/CDUC_species.tsv" -y "$input_dir/CDUC_EggNOG.tsv" \
    -o "$output_root/halla_EggNOG_CDUC" \
    -m spearman --fdr_alpha 0.05 --permute_iters 1000

halla -x "$input_dir/CDUC_species.tsv" -y "$input_dir/CDUC_pathway.tsv" \
    -o "$output_root/halla_pathway_CDUC" \
    -m spearman --fdr_alpha 0.05 --permute_iters 1000
