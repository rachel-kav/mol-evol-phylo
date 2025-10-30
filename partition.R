#partition sites in segments

library(ape)
library(Biostrings)
library(seqinr)

fasta <- "alignments/trim_genome_50.fasta"
aln_mat <- read.dna(fasta, format = "fasta") 
aln_length <- ncol(aln_mat)
num_seq <- nrow(aln_mat)
window_size <- round(aln_length / 2500) #pick number of segments
n_windows <- floor(aln_length / window_size)
window_size <- aln_length / n_windows

trees <- list()

dir.create("iqtree_output_50", showWarnings = FALSE)

for (start in seq(1, aln_length, by = window_size)) {
  end <- start + window_size - 1
  print(start)
  
  if (end > aln_length) {
    end <- aln_length
  } #Make sure the last segment is not out of bounds
  print(end)
  #current window 
  segment <- aln_mat[, start:end]
  
  file_path <- paste0("iqtree_output_50/aln_segment_", start, "_to_", end)
  write.dna(segment, file = paste0(file_path, ".fasta"), format = "fasta", colsep = "")


  system(paste0("iqtree2 -s ", file_path, ".fasta -m GTR -nt AUTO -fast -pre ", file_path))
  
  treefile <- paste0(file_path, ".treefile")
  if (file.exists(treefile)) {
    trees[[paste0("Tree_", start, "_to_", end)]] <- read.tree(treefile)
  } else {
    cat("Skipping: no treefile for", start, "to", end, "\n")
  }
}

if (length(trees) > 0) {
  write.tree(do.call(c, trees), file = "iqtree_output_50/trim_trees.newick")
  cat("Wrote combined tree file: trim_trees.newick\n")
} else {
  cat("No trees generated.\n")
}

