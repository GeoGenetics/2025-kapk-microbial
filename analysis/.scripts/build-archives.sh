#!/usr/bin/env bash
# Build data.tar.gz and results.tar.gz for ERDA distribution.
# Run from analysis/ directory.
set -euo pipefail
cd "$(dirname "$0")/.."

DATA_SRC=$(readlink -f data)  # resolve symlink

echo "=== Building data.tar.gz ==="
# All data files except:
#   imgvr_gene_ranks.tsv       — 272MB pipeline intermediate, not needed for R scripts
#   methanoflorens/            — 30GB; we include only the specific files needed by script 13
#   viruses-aa IMGVR4 files    — permission denied, not used by any current script

# Step 1: main data (excluding methanoflorens entirely)
tar czf data.tar.gz \
    --exclude='./methanoflorens' \
    --exclude='./cdata/imgvr_gene_ranks.tsv' \
    --exclude='./taxonomy/viruses-aa/IMGVR4_derepG-archaea-profiling.group-abundances-agg.tsv.gz' \
    --exclude='./taxonomy/viruses-aa/IMGVR4_derepG-archaea-profiling.group-abundances.tsv.gz' \
    -C "$DATA_SRC" \
    --transform 's|^\./|data/|' \
    .

# Step 2: build methanoflorens file list and append to archive
META="$DATA_SRC/methanoflorens"
METHANO_FILES=(
    beast2_nodec/tipdating/codon_tipdating_mcc.tree
    beast2_nodec/tipdating/codon_tipdating.log
    beast2_nodec/2clock_stasis/codon_2clock_stasis_nodec.log
    beast2_nodec/ucld_tipdating/codon_tipdating.log
    beast2_nodec/freedate_rep3/codon_freedate_nodec.log
    beast2_nodec/freedate_rep4/codon_freedate_nodec.log
    phylogenomics_79taxa_bak/08_trees/codon12_ml_tree.treefile
    phylogenomics_79taxa_bak/11_focal_snps/focal_pairwise_snp_counts.tsv
    phylogenomics_79taxa_bak/11_focal_snps/focal_pairwise_ani.tsv
    fastani/deconvolved_vs_refs.tsv
    deconvolve/read_length_histogram.tsv
    deconvolve/damage_model.tsv
    deconvolve/smiley_data.tsv
    s17_assembly/deconvolve/read_length_histogram.tsv
    s17_assembly/deconvolve/smiley_data.tsv
)
# Add ps_*/beast_ps.out files
while IFS= read -r f; do
    rel="${f#$META/}"
    METHANO_FILES+=("$rel")
done < <(find "$META/beast2_nodec" -name "beast_ps.out" 2>/dev/null)

echo "Merging methanoflorens subset into data.tar.gz..."
mkdir -p /tmp/data_merge/data/methanoflorens
for f in "${METHANO_FILES[@]}"; do
    dest="/tmp/data_merge/data/methanoflorens/$f"
    mkdir -p "$(dirname "$dest")"
    cp "$META/$f" "$dest"
done
# Repack: extract existing, add methanoflorens, repack
tar xzf data.tar.gz -C /tmp/data_merge
tar czf data.tar.gz -C /tmp/data_merge data/
rm -rf /tmp/data_merge
echo "data.tar.gz: $(du -sh data.tar.gz | cut -f1)"

echo ""
echo "=== Building results.tar.gz ==="
RES=results

{
    find $RES/controls $RES/taxonomy $RES/damage $RES/sourcetracker $RES/figures \
        -type f 2>/dev/null
    find $RES/functional_agp -type f \( \
        -name "anvio_modules.txt" -o \
        -name "anvio_hits.txt" -o \
        -name "emi.functional.tsv" -o \
        -name "emi.protein.slim.tsv.gz" -o \
        -name "pathway_damage_stats.tsv" -o \
        -name "cazy_emi.functional.tsv" -o \
        -name "cazy_damage.slim.tsv.gz" -o \
        -name "kegg_module_damage.tsv" -o \
        -name "cazy_family_damage.tsv" \
    \) 2>/dev/null
    find $RES/virome_agp/functional -name "viral_emi.functional.tsv" 2>/dev/null
    echo "$RES/virome_agp/viral_auth.per_sample.tsv"
} > /tmp/results_files.txt

tar czf results.tar.gz --files-from=/tmp/results_files.txt
rm /tmp/results_files.txt
echo "results.tar.gz: $(du -sh results.tar.gz | cut -f1)"
echo "Done."
