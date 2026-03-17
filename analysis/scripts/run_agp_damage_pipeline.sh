#!/bin/bash
# run_agp_damage_pipeline.sh
#
# Runs the AGP damage annotation pipeline for functional profiling.
# This generates .protein_damage.tsv files that are consumed by 07a--agp-functional.R
#
# Prerequisites:
#   - AGP installed and in PATH (or specify AGP_BIN)
#   - MMseqs2 installed
#   - VTML20 substitution matrix available
#
# Usage:
#   ./run_agp_damage_pipeline.sh <sample_id> <input_fastq> <db_name>
#
# Example:
#   ./run_agp_damage_pipeline.sh SampleA data/reads/SampleA.fq.gz kegg

set -euo pipefail

# Configuration
AGP_BIN="${AGP_BIN:-/vol/cloud/agp/bin/agp}"
MMSEQS="${MMSEQS:-mmseqs}"
VTML20_MATRIX="${VTML20_MATRIX:-/maps/projects/fernandezguerra/apps/repos/read2Struct/read2Struct/data/mmseqs2/sub-mat/VTML20.out}"
THREADS="${THREADS:-8}"

# Database paths (adjust as needed)
KEGG_DB="/datasets/fernandezguerra/data/DB/KEGG/20220228/genes/mmseqs2/kegg_genes-i0.9-c0.8"
CAZY_DB="${CAZY_DB:-./data/db/cazy_mmseqs}"
VIRAL_DB="${VIRAL_DB:-./data/db/viral_mmseqs}"

# Parse arguments
if [[ $# -lt 3 ]]; then
    echo "Usage: $0 <sample_id> <input_fastq> <db_name>"
    echo ""
    echo "Arguments:"
    echo "  sample_id    - Sample identifier"
    echo "  input_fastq  - Path to input FASTQ file (can be gzipped)"
    echo "  db_name      - Database to search: kegg, cazy, or viral"
    echo ""
    echo "Environment variables:"
    echo "  AGP_BIN       - Path to AGP binary (default: /vol/cloud/agp/bin/agp)"
    echo "  MMSEQS        - Path to mmseqs binary (default: mmseqs)"
    echo "  VTML20_MATRIX - Path to VTML20 substitution matrix"
    echo "  THREADS       - Number of threads (default: 8)"
    exit 1
fi

SAMPLE_ID="$1"
INPUT_FASTQ="$2"
DB_NAME="$3"

# Select database
case "$DB_NAME" in
    kegg)
        TARGET_DB="$KEGG_DB"
        ;;
    cazy)
        TARGET_DB="$CAZY_DB"
        ;;
    viral)
        TARGET_DB="$VIRAL_DB"
        ;;
    *)
        echo "Error: Unknown database '$DB_NAME'. Use: kegg, cazy, or viral"
        exit 1
        ;;
esac

# Output directories
OUT_DIR="./results/agp_damage/${SAMPLE_ID}/${DB_NAME}"
TMP_DIR="${OUT_DIR}/tmp"
mkdir -p "$OUT_DIR" "$TMP_DIR"

echo "=== AGP Damage Annotation Pipeline ==="
echo "Sample:   $SAMPLE_ID"
echo "Input:    $INPUT_FASTQ"
echo "Database: $DB_NAME ($TARGET_DB)"
echo "Output:   $OUT_DIR"
echo ""

# Step 1: AGP predict with damage index
echo "[1/3] Running AGP predict..."
GFF_OUT="${OUT_DIR}/${SAMPLE_ID}.gff"
FAA_OUT="${OUT_DIR}/${SAMPLE_ID}.search.faa"
AGD_OUT="${OUT_DIR}/${SAMPLE_ID}.agd"

if [[ ! -f "$AGD_OUT" ]]; then
    "$AGP_BIN" predict \
        -i "$INPUT_FASTQ" \
        -o "$GFF_OUT" \
        --fasta-aa-masked "$FAA_OUT" \
        --damage-index "$AGD_OUT" \
        --adaptive \
        -t "$THREADS" \
        2>&1 | tail -20
    echo "  Created: $AGD_OUT"
else
    echo "  Skipping (exists): $AGD_OUT"
fi

# Step 2: MMseqs2 search with VTML20
echo ""
echo "[2/3] Running MMseqs2 search..."
HITS_OUT="${OUT_DIR}/${SAMPLE_ID}_${DB_NAME}_hits.tsv"

if [[ ! -f "$HITS_OUT" ]]; then
    "$MMSEQS" easy-search \
        "$FAA_OUT" \
        "$TARGET_DB" \
        "$HITS_OUT" \
        "$TMP_DIR" \
        --search-type 1 \
        --min-length 12 \
        -e 10.0 \
        --min-seq-id 0.86 \
        -c 0.65 \
        --cov-mode 2 \
        --sub-mat "$VTML20_MATRIX" \
        --seed-sub-mat "$VTML20_MATRIX" \
        -s 2 \
        -k 6 \
        --spaced-kmer-pattern 11011101 \
        --threads "$THREADS" \
        --format-output "query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qlen,tlen,qaln,taln" \
        2>&1 | tail -10
    echo "  Created: $HITS_OUT"
else
    echo "  Skipping (exists): $HITS_OUT"
fi

# Step 3: AGP damage-annotate
echo ""
echo "[3/3] Running AGP damage-annotate..."
DAMAGE_OUT="${OUT_DIR}/${SAMPLE_ID}.protein_damage.tsv"

if [[ ! -f "$DAMAGE_OUT" ]]; then
    "$AGP_BIN" damage-annotate \
        -i "$HITS_OUT" \
        --damage-index "$AGD_OUT" \
        -o "$DAMAGE_OUT" \
        2>&1 | tail -10
    echo "  Created: $DAMAGE_OUT"
else
    echo "  Skipping (exists): $DAMAGE_OUT"
fi

# Cleanup temp files
rm -rf "$TMP_DIR"

echo ""
echo "=== Pipeline Complete ==="
echo "Damage annotations: $DAMAGE_OUT"
echo ""
echo "Next steps:"
echo "  1. Run for all samples"
echo "  2. Execute: Rscript scripts/07a--agp-functional.R"
