
# set your working directory to where the .cf.stat files are
setwd("/Users/morah/Library/CloudStorage/OneDrive-ITU/Documents/Github repos/Courses/Molecular Evolution and Phylogenomic course/mol-evol-phylo")

# list all .cf.stat files
files <- list.files(path = "iqtree_output_africa_scf", pattern = "\\.cf\\.stat$", full.names = TRUE)

# read each file and compute mean sCF
results <- data.frame(
  file = files,
  mean_sCF = sapply(files, function(f) {
    x <- read.table(f, header = TRUE, sep = "\t", fill = TRUE, stringsAsFactors = FALSE)
    mean(x$sCF, na.rm = TRUE)
  })
)

# print the results
print(results)

# load the MAST weights summary file
mast_weights <- read.csv("mast_africa_output/mast_weights_summary_sort.csv")

# Extract segment names from the file paths in results
results$segment_name <- gsub(".*/(.*?)_scf\\.cf\\.stat$", "\\1", results$file)

# Extract the starting position from each interval and sort
results$start_pos <- as.numeric(gsub(".*_(\\d+)_to_\\d+", "\\1", results$segment_name))
results <- results[order(results$start_pos), ]

# Assign sequential tree numbers based on sorted order
results$tree_nr <- paste0("tree_", seq_len(nrow(results)))

# Remove .treefile extension from mast_weights tree names for matching
mast_weights$tree_name_clean <- gsub("\\.treefile$", "", mast_weights$tree_name)

# Merge the dataframes based on matching segment names
mast_weights_with_scf <- merge(mast_weights, results[, c("segment_name", "mean_sCF", "tree_nr")], 
                               by.x = "tree_name_clean", by.y = "segment_name", 
                               all.x = TRUE)

# Display the merged dataframe
print(mast_weights_with_scf)

# Sort by weight in descending order
mast_weights_with_scf <- mast_weights_with_scf[order(mast_weights_with_scf$weight, decreasing = TRUE), ]

# Remove the aln_segment_142409_to_142412 entry
mast_weights_with_scf <- mast_weights_with_scf[mast_weights_with_scf$tree_name_clean != "aln_segment_142409_to_142412", ]

# Save plot to PNG
png("MAST_weights_and_sCF_plot.png", width = 1200, height = 800, res = 100)

# Increase bottom margin to accommodate rotated labels
par(mar = c(18, 4, 4, 4))

# Create the barplot of weights (left y-axis)
bp <- barplot(mast_weights_with_scf$weight, 
              names.arg = mast_weights_with_scf$tree_nr, 
              main = "MAST Weights and sCF Values", 
              ylab = "MAST Weight",
              las = 2,
              col = "lightblue")

# Add a second y-axis for sCF values
par(new = TRUE)

# Plot sCF values as points/lines (right y-axis)
plot(bp, mast_weights_with_scf$mean_sCF, 
     type = "p", 
     pch = 19, 
     col = "red", 
     cex = 1.2,
     axes = FALSE, 
     xlab = "", 
     ylab = "",
     ylim = c(40, max(mast_weights_with_scf$mean_sCF) * 1.1))

# Add the right y-axis for sCF
axis(4, col = "red", col.axis = "red")
mtext("Mean sCF", side = 4, line = 2.5, col = "red")

# Add lines connecting the sCF points
lines(bp, mast_weights_with_scf$mean_sCF, col = "red", lwd = 2)

# Add a legend
legend("topright", 
       legend = c("MAST Weight", "Mean sCF"), 
       fill = c("lightblue", NA),
       pch = c(NA, 19),
       col = c("lightblue", "red"),
       bty = "n")

# Close and save the PNG
dev.off()

# Reset margins to default
par(mar = c(5, 4, 4, 2))

# Print confirmation message
cat("Plot saved as 'MAST_weights_and_sCF_plot.png' in the current working directory\n")


