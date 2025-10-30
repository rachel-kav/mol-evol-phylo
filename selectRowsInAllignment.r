if (!require(ape)) install.packages("ape", repos = "https://cloud.r-project.org")
library(ape)

# ---- User input ----
# Path to input and output FASTA files
input_fasta  <- "your_sequences.fasta"          # change this
output_fasta <- "filtered_sequences.fasta"      # change this

# List of taxa (sequence names) to remove
remove_taxa <- c("Outgroup", "HIV2_Ref")        # edit as needed

# Optionally, remove sequences that match a pattern (e.g. any name containing "HIV2")
pattern_to_remove <- "HIV2"                     # set to "" to disable

# ---- Read FASTA ----
cat("Reading FASTA file:", input_fasta, "\n")
fasta <- read.FASTA(input_fasta)
cat("Number of sequences read:", length(fasta), "\n")

# ---- Inspect names ----
cat("Example sequence names:\n")
print(head(names(fasta)))

# ---- Remove exact matches ----
if (length(remove_taxa) > 0) {
  fasta <- fasta[!names(fasta) %in% remove_taxa]
  cat("Removed taxa:", paste(remove_taxa, collapse = ", "), "\n")
}

# ---- Remove by pattern ----
if (pattern_to_remove != "") {
  fasta <- fasta[!grepl(pattern_to_remove, names(fasta))]
  cat("Removed taxa matching pattern:", pattern_to_remove, "\n")
}

# ---- Write filtered FASTA ----
write.FASTA(fasta, file = output_fasta)
cat("Filtered FASTA written to:", output_fasta, "\n")
cat("Number of sequences remaining:", length(fasta), "\n")