#!/bin/bash

# List of FASTA files
FASTAS=("alignments/UL15_non_recomb.fst" "alignments/UL15_recomb.fst" "alignments/UL29_non_recomb.fst" "alignments/UL29_recomb.fst" "alignments/UL30_non_recomb.fst" "alignments/UL30_recomb.fst" "alignments/UL39_non_recomb.fst" "alignments/UL39_recomb.fst")

# Make output folder
mkdir -p iqtree_genes

for FASTA in "${FASTAS[@]}"; do
    # Extract base name without folder and extension
    BASENAME=$(basename "$FASTA" .fst)
    OUT1="iqtree_genes/$BASENAME"

    echo "Running IQ-TREE2 for $FASTA ..."
    iqtree2 -s "$FASTA" -m "GTR+FO+R4" --prefix "$OUT1" -fast -nt AUTO -redo
done
