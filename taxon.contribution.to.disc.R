# Requires package TreeDist - Takes the distance between trees where each taxon has been removed as the discordance metric for the taxon.
# The greater the distance, the greater the contribution of a taxon to tree discordance.

library(TreeDist)
library(ape)
library(ggplot2)

taxon.contribution.to.disc <- function(tr1, tr2){
	taxa <- intersect(tr1$tip.label, tr2$tip.label)
	disc.metric <- vector()
	for(i in 1:length(taxa)){
		# For an inverse, where higher values are taxa that contribute least to discordance
		#disc.metric[i] <- ClusteringInfoDistance(drop.tip(tr1, taxa[i]), drop.tip(tr2, taxa[i]))
		# More intuitive, where higher values are taxa that contribute most to discordance 
		disc.metric[i] <- MutualClusteringInfo(drop.tip(tr1, taxa[i]), drop.tip(tr2, taxa[i]))
	}
	names(disc.metric) <- taxa
	return(disc.metric)
}

# Load all .treefile files from the iqtree_tree_africa directory
tree_files <- list.files("iqtree_tree_africa", pattern = "\\.treefile$", full.names = TRUE)
trees <- lapply(tree_files, read.tree)
names(trees) <- basename(tree_files)
# Visualize the first tree
plot(trees[[1]], main = names(trees)[1])

# Example usage of taxon.contribution.to.disc function
tr1 <- trees[[1]]
tr2 <- trees[[8]]


plot(tr1, main = names(tr1))
plot(tr2, main = names(tr2))
disc_metrics <- taxon.contribution.to.disc(tr1, tr2)
print(disc_metrics)
sort(disc_metrics, decreasing = TRUE)
barplot(sort(disc_metrics, decreasing = TRUE), las=2, main="Taxon Contribution to Discordance between Tree 1 and Tree 2", ylab="Discordance Metric", ylim=c(min(disc_metrics)-0.5, max(disc_metrics)*1.1))



### Forget about this
# Create pairwise comparisons for all tree combinations (upper triangle + diagonal)
n_trees <- length(trees)

# Save to PDF file
pdf("tree_comparisons.pdf", width = 16, height = 14)

par(mfrow = c(n_trees, n_trees), mar = c(1, 1, 2, 0.5))

for(i in 1:n_trees) {
	for(j in 1:n_trees) {
		if(i == j) {
			# Diagonal: show the tree itself
			plot(trees[[i]], main = paste("Tree", i))
		} else if (i < j) {
			# Upper triangle: show discordance barplot
			disc_metrics <- taxon.contribution.to.disc(trees[[i]], trees[[j]])
			barplot(sort(disc_metrics, decreasing = TRUE), 
							main = paste("Tree", i, "vs Tree", j), 
							names.arg = "",
							las = 2, cex.axis = 0.6, ylim = c(min(disc_metrics)*0.8, max(disc_metrics) * 1.1))
		} else {
			# Lower triangle: leave empty (plot nothing)
			plot.new()
		}
	}
}

# Reset plotting parameters
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2))

# Close and save the PDF file
dev.off()

# Print message about where file was saved
cat("Plot saved as 'tree_comparisons.pdf' in the current working directory\n")

### Forget about this part above



### This analysis is done for initial trees
# Get taxon names from the first tree (assuming all trees have the same taxa)
taxon_names <- trees[[1]]$tip.label
taxonContribution <- array(0, dim = c(length(taxon_names)))
names(taxonContribution) <- taxon_names

for(i in 1:n_trees) {
	for(j in 1:n_trees) 
		if (i < j) {
		{
			# Upper triangle: show discordance barplot
			disc_metrics <- taxon.contribution.to.disc(trees[[i]], trees[[j]])
			taxonContribution <- taxonContribution + disc_metrics
		} 
	}
}
taxonContribution <- taxonContribution / (n_trees*(n_trees-1)/2)


sort(taxonContribution, decreasing = TRUE)

# Create color vector based on HIV type
colors <- ifelse(grepl("Human_herpesvirus_1", names(sort(taxonContribution, decreasing = TRUE))), "blue", "red")

barplot(sort(taxonContribution, decreasing = TRUE), 
				main = "Average Taxon Contribution to Discordance across all Tree Pairs", 
				las = 2, cex.axis = 0.6, ylim = c(min(taxonContribution)*0.8, max(taxonContribution) * 1.1),
				col = colors)

# Add legend
legend("topright", legend = c("Human_herpesvirus_1", "Human_herpesvirus_2"), fill = c("blue", "red"))

### this analysis is done for mast optimised trees
# Load MAST weights summary
mast_weights <- read.csv("mast_africa_output/mast_weights_summary_sort.csv")
# Load MAST optimized tree
mast_trees <- read.tree("mast_africa_output/mast1.treefile")

# Assign specific segment names to MAST trees
segment_names <- c("aln_segment_1_to_10172",
                   "aln_segment_10173_to_20344",
                   "aln_segment_20345_to_30516",
                   "aln_segment_30517_to_40688",
                   "aln_segment_40689_to_50860",
                   "aln_segment_50861_to_61032",
                   "aln_segment_61033_to_71204",
                   "aln_segment_71205_to_81376",
                   "aln_segment_81377_to_91548",
                   "aln_segment_91549_to_101720",
                   "aln_segment_101721_to_111892",
                   "aln_segment_111893_to_122064",
                   "aln_segment_122065_to_132236",
                   "aln_segment_132237_to_142408",
                   "aln_segment_142409_to_142412")

names(mast_trees) <- segment_names
# Access individual trees from the multitree file
# Order mast_trees according to weights in mast_weights
# Remove .treefile extension from mast_weights names for matching
mast_weights$tree_name_clean <- gsub("\\.treefile$", "", mast_weights$tree_name)

# Reorder mast_trees according to this ordering
ordered_trees <- list()
for (j in 1:length(mast_weights$tree_name_clean)) {
	ordered_trees[[j]] <- mast_trees[[mast_weights$tree_name_clean[j]]]
}
names(ordered_trees) <- mast_weights$tree_name_clean

for (j in 2:3) {
	taxon_names <- ordered_trees[[j]]$tip.label
	taxonContribution <- array(0, dim = c(length(taxon_names)))
	names(taxonContribution) <- taxon_names
	for(i in 1:j) {
		for(k in 1:j) {
			if (i < k) {
				{
				# Upper triangle: show discordance barplot
				disc_metrics <- taxon.contribution.to.disc(ordered_trees[[i]], ordered_trees[[k]])
				taxonContribution <- taxonContribution + disc_metrics
				} 
			}
		}
	}
taxonContribution <- taxonContribution / (j*(j-1)/2)
sort(taxonContribution, decreasing = TRUE)

# Create color vector based on HIV type
colors <- ifelse(grepl("Human_herpesvirus_1", names(sort(taxonContribution, decreasing = TRUE))), "blue", "red")

# Save the plot to PDF in the subfolder
pdf(paste0("discordanceContributionPlots/mast_taxon_contribution_discordance_w_", j, "trees.pdf"), width = 10, height = 8)

barplot(sort(taxonContribution, decreasing = TRUE), 
				main = paste("Average Taxon Contribution to Discordance across", j, "best Tree Pairs"), 
				las = 2, cex.axis = 0.6, ylim = c(min(taxonContribution)*0.8, max(taxonContribution) * 1.1),
				col = colors)

legend("topright", legend = c("Human_herpesvirus_1", "Human_herpesvirus_2"), fill = c("blue", "red"))

dev.off()

cat("Plot saved as 'discordanceContributionPlots/mast_taxon_contribution_discordance_w_", j, "trees.pdf'\n")
}
# Create a multiplot showing taxon contribution to discordance for three tree combinations
pdf("discordanceContributionPlots/three_tree_combinations.pdf", width = 15, height = 5)

par(mfrow = c(1, 3), mar = c(8, 4, 4, 2))

# Define the three combinations
combinations <- list(c(1, 2), c(1, 3), c(2, 3))

# Create consistent color palettes for each virus type
hiv1_taxa <- names(sort(disc_metrics, decreasing = TRUE))[grepl("Human_herpesvirus_1", names(sort(disc_metrics, decreasing = TRUE)))]
hiv2_taxa <- names(sort(disc_metrics, decreasing = TRUE))[grepl("Human_herpesvirus_2", names(sort(disc_metrics, decreasing = TRUE)))]

# Generate color palettes
blue_palette <- colorRampPalette(c("#000080", "#87CEEB"))(length(hiv1_taxa))
red_palette <- colorRampPalette(c("#8B0000", "#FFB6C1"))(length(hiv2_taxa))

# Create named color vectors for consistency
hiv1_colors <- setNames(blue_palette, hiv1_taxa)
hiv2_colors <- setNames(red_palette, hiv2_taxa)

for (combo_idx in 1:length(combinations)) {
	i <- combinations[[combo_idx]][1]
	j <- combinations[[combo_idx]][2]
	
	# Calculate discordance for this combination
	disc_metrics <- taxon.contribution.to.disc(ordered_trees[[i]], ordered_trees[[j]])
	sorted_names <- names(sort(disc_metrics, decreasing = TRUE))
	
	# Create color vector with consistent nuances
	colors <- character(length(sorted_names))
	for (k in 1:length(sorted_names)) {
		if (sorted_names[k] %in% names(hiv1_colors)) {
			colors[k] <- hiv1_colors[sorted_names[k]]
		} else if (sorted_names[k] %in% names(hiv2_colors)) {
			colors[k] <- hiv2_colors[sorted_names[k]]
		} else {
			colors[k] <- "gray"  # fallback color
		}
	}
	
	barplot(sort(disc_metrics, decreasing = TRUE), 
					main = paste("RF between ", i, " most contributing tree and ", j," most contributing tree", sep=""), 
					las = 2, cex.axis = 0.6, ylim = c(min(disc_metrics), max(disc_metrics) * 1.05),
					col = colors, cex.names = 0.5, names.arg = sorted_names)
	
	if (combo_idx == 1) {
		legend("topright", legend = c("Human_herpesvirus_1", "Human_herpesvirus_2"), 
					 fill = c("blue", "red"), cex = 0.8)
	}
}

par(mfrow = c(1, 1), mar = c(5, 4, 4, 2))
dev.off()

cat("Multiplot saved as 'discordanceContributionPlots/three_tree_combinations.pdf'\n")


###





