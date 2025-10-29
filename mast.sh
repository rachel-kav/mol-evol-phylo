#!/bin/bash

#-s sequence alignment file
#-m substitution model
#+T adds MAST
#-te tree output
FASTA="alignments/trim_genome_50.fasta"
TREES="iqtree_output_50/trim_trees.newick"
OUT1="mast_output/trim_genome_50"

iqtree2 -s $FASTA -m "GTR+FO+R4+T" -te $TREES --prefix $OUT1 -fast -nt AUTO -wspmr -wslmr -redo