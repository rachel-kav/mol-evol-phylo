# Repeat MAST iteratively by tree weight
library(ape)
library(future.apply)

setwd("/home/rachel/mol_evol")

# === 1. Load tree weights from MAST ===
weights_df <- read.csv("iqtree_output_africa/mast_weights_summary_sort.csv")

# === 2. Load and reorder trees ===
nwk <- read.tree("iqtree_output_africa/trim_trees.newick")
tree_files <- basename(weights_df$tree_name)

# === 3. Create newick subsets ordered by tree weight ===
dir.create("tree_subsets", showWarnings = FALSE, recursive = TRUE)
for (i in 2:nrow(weights_df)) {
  subset_trees <- nwk[match(weights_df$tree_name[1:i], tree_files)]
  write.tree(subset_trees, file.path("tree_subsets", paste0("modeltrees_final_", i, ".newick")))
}

# === 4. Define paths ===
FASTA <- "trim_genome_africa.fasta"
TREES_SUBSETS <- "tree_subsets"
OUT <- "output_mast_repeats"

# === 5. Build IQ-TREE commands ===
commands <- sapply(2:nrow(weights_df), function(i) {
  paste0(
    "mkdir -p ", OUT, "/mast_model_", i, " && ",
    "iqtree2 -s ", FASTA,
    " -m 'GTR+FO+R2+T'",
    " -te ", TREES_SUBSETS, "/modeltrees_final_", i, ".newick",
    " --prefix ", OUT, "/mast_model_", i, "/mast_model_", i,
    " -nt 2 -wslmr -redo"
  )
})

# === 6. Run in parallel ===
plan(multisession, workers = 6)
future_lapply(commands, function(cmd) {
  message("Running: ", cmd)
  system(cmd)
}, future.seed = TRUE)
