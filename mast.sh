#!/bin/bash

#-s sequence alignment file
#-m substitution model
#+T adds MAST
#-te tree output
#cd /home/rachel/mol_evol

FASTA="trim_genome_50.fasta"
TREES="iqtree_output_50/trim_trees.newick"
OUT1="mast_output/mast1"

iqtree2 -s $FASTA -m "GTR+FO+R4+T" -te $TREES --prefix $OUT1 -fast -nt AUTO -wslmr -redo
#wslmr calculates the site log marginal likelihoods