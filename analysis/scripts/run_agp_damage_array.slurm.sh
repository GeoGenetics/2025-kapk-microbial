#!/bin/bash
#SBATCH --job-name=agp_damage
#SBATCH --output=agp_damage_%A_%a.out
#SBATCH --error=agp_damage_%A_%a.err
#SBATCH --time=4:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8
#SBATCH --array=0-49%10

# run_agp_damage_array.slurm.sh
#
# SLURM array job for running AGP damage annotation on all KapK samples.
#
# Usage:
#   # First, generate sample list:
#   ls data/reads/*.fq.gz | xargs -I{} basename {} .fq.gz > sample_list.txt
#
#   # Adjust array range based on number of samples:
#   wc -l sample_list.txt  # e.g., 50 samples -> --array=0-49
#
#   # Submit:
#   sbatch scripts/run_agp_damage_array.slurm.sh

set -euo pipefail

# Configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORK_DIR="$(dirname "$SCRIPT_DIR")"
cd "$WORK_DIR"

# Sample list file
SAMPLE_LIST="${SAMPLE_LIST:-./sample_list.txt}"
READS_DIR="${READS_DIR:-./data/reads}"
DB_NAME="${DB_NAME:-kegg}"

# Get sample for this array task
SAMPLE_ID=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$SAMPLE_LIST")

if [[ -z "$SAMPLE_ID" ]]; then
    echo "Error: No sample found for task $SLURM_ARRAY_TASK_ID"
    exit 1
fi

# Input file
INPUT_FASTQ="${READS_DIR}/${SAMPLE_ID}.fq.gz"

if [[ ! -f "$INPUT_FASTQ" ]]; then
    # Try alternative extensions
    for ext in ".fastq.gz" ".fq" ".fastq"; do
        alt="${READS_DIR}/${SAMPLE_ID}${ext}"
        if [[ -f "$alt" ]]; then
            INPUT_FASTQ="$alt"
            break
        fi
    done
fi

if [[ ! -f "$INPUT_FASTQ" ]]; then
    echo "Error: Input file not found for sample $SAMPLE_ID"
    echo "Tried: ${READS_DIR}/${SAMPLE_ID}.{fq.gz,fastq.gz,fq,fastq}"
    exit 1
fi

echo "=== SLURM Array Task $SLURM_ARRAY_TASK_ID ==="
echo "Sample: $SAMPLE_ID"
echo "Input:  $INPUT_FASTQ"
echo "DB:     $DB_NAME"
echo ""

# Set environment
export AGP_BIN="${AGP_BIN:-/vol/cloud/agp/bin/agp}"
export THREADS=$SLURM_CPUS_PER_TASK

# Run pipeline
./scripts/run_agp_damage_pipeline.sh "$SAMPLE_ID" "$INPUT_FASTQ" "$DB_NAME"

echo ""
echo "Task $SLURM_ARRAY_TASK_ID complete."
