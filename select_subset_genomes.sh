#!/bin/bash
set -euo pipefail

INPUT_DIR="genomes/ncbi_dataset/data"
OUTPUT_DIR="subset_genomes"
NUM_GENOMES=10

mkdir -p "$OUTPUT_DIR"

echo "[Step 1] Listing all genome folders..."
GENOME_FOLDERS=($(find "$INPUT_DIR" -maxdepth 1 -type d -name "GCF_*"))

TOTAL=${#GENOME_FOLDERS[@]}
echo "Found $TOTAL genome folders."

if [[ $TOTAL -lt $NUM_GENOMES ]]; then
    echo "❌ Not enough genomes to select $NUM_GENOMES. Found only $TOTAL."
    exit 1
fi

echo "[Step 2] Randomly selecting $NUM_GENOMES genomes..."
# Shuffle and pick first $NUM_GENOMES
SHUFFLED=($(printf "%s\n" "${GENOME_FOLDERS[@]}" | shuf | head -n $NUM_GENOMES))

echo "[Step 3] Copying selected genomes to $OUTPUT_DIR..."
for GENOME in "${SHUFFLED[@]}"; do
    BASENAME=$(basename "$GENOME")
    echo "Copying $BASENAME..."
    cp -r "$GENOME" "$OUTPUT_DIR/"
done

echo "[Step 4] Done! $NUM_GENOMES genomes copied to $OUTPUT_DIR/"
