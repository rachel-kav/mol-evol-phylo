if (!require(ape)) install.packages("ape", repos = "https://cloud.r-project.org")
library(ape)
setwd('/Users/morah/Library/CloudStorage/OneDrive-ITU/Documents/Github repos/Courses/Molecular Evolution and Phylogenomic course/mol-evol-phylo/')
# ---- User input ----
# Path to input and output FASTA files
input_fasta  <- "alignments/Whole_Genome_1.fst"          # change this
output_fasta <- "filtered_sequences.fasta"      # change this

# ---- Read FASTA ----
cat("Reading FASTA file:", input_fasta, "\n")
fasta <- read.FASTA(input_fasta)
cat("Number of sequences read:", length(fasta), "\n")

# ---- Inspect names ----
cat("Example sequence names:\n")
print(head(names(fasta)))

# ---- Remove exact matches ----
keep_patterns <- c(
  "HM585496",
  "HM585497",
  "HM585498",
  "HM585499",
  "HM585500",
  "HM585501",
  "HM585502",
  "HM585503",
  "HM585504",
  "HM585505",
  "HM585506",
  "HM585507",
  "HM585509",
  "HM585510",
  "HM585511",
  "KR135299",
  "KR135300",
  "KR135301",
  "KR135302",
  "KR135303",
  "KR135304",
  "KR135305",
  "KR135306",
  "KR135307",
  "KR135315",
  "KR135316",
  "KR135317",
  "KR135318",
  "KR135319"
)

# Combine into one regular expression pattern
pattern <- paste(keep_patterns, collapse = "|")

# Logical vector: TRUE if name contains any of those country names
keep_idx <- grepl(pattern, names(fasta), ignore.case = TRUE)

# Subset
filtered_fasta <- fasta[keep_idx]


# ---- Write filtered FASTA ----
write.FASTA(filtered_fasta, file = output_fasta)
cat("Filtered FASTA written to:", output_fasta, "\n")


