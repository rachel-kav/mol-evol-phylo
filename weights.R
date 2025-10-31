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

mast_summary_sort <- read.csv("mast_africa_output/mast_weights_summary_sort.csv")

#BIC values

data <- read.table("output_mast_repeats/bic_summary.tsv", header = TRUE, sep = "\t")

#data_sorted <- data %>%
#  arrange(desc(bic))

data$model_num <- as.numeric(sub("mast_model_", "", data$model))
data_sorted <- data[order(data$model_num), ]

BIC_table <- data_sorted %>% mutate("BIC_transformed" = bic/10000)

bic_plot <- ggplot(BIC_table, aes(x = factor(model_num), y = BIC_transformed)) +
  geom_point(colour = "black", size = 3) +
  theme_minimal() +
  labs(x = "Number of Trees in MAST Model",
       y = "BIC (x10e5)")
#order and add BIC curve

ggsave("figures/bic_plot.png", bic_plot, width = 10, height = 6)

ggsave("figures/mast_weights_plot.png", plot, width = 10, height = 6)


