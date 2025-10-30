library(ape)
library(Biostrings)
library(tidyverse)
library(ggplot2)

setwd("/Users/qjs599/Desktop/mol-evol-phylo")
trees <- read.tree("iqtree_tree_africa/trim_trees.newick")

tree_files <- list.files("iqtree_tree_africa/", pattern = "*.treefile$", full.names = TRUE)

tree_names <- basename(tree_files)

weights <- c(0.04653,0.04049,0.03409,0.17885,0.12364,
             0.26024,0.03934,0.04039,0.00844,0.07515,
             0.06122,0.02150,0.02636,0.04363,0.00015)

mast_summary <- data.frame(
  tree_name = tree_names,
  weight = weights
)

mast_summary <- mast_summary %>%
  mutate(start = as.numeric(sub("aln_segment_([0-9]+)_.*", "\\1", tree_name))) %>%
  arrange(start) %>%
  select(-start)

mast_summary_sort <- mast_summary %>%
  arrange(desc(weight))

write.csv(mast_summary_sort, "mast_africa_output/mast_weights_summary_sort.csv", row.names = FALSE)

p <- ggplot(mast_summary_sort, aes(x = tree_name, y = weight)) +
        geom_col(position = "dodge") +
        xlab("Tree") +
        ylab("Weight") +
        theme_bw()

ggplot(mast_summary, aes(x = tree_name, y = weight)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))