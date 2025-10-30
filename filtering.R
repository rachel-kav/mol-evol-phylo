library(Biostrings)
library(ape)

trimCols <- function(al, prop, codon = T){
	 mat <- as.character(as.matrix(al))
	 ntax <- nrow(mat)
	 propthres <- 1-prop
	 compliantSites <- apply(mat, 2, function(x){
	 		x <- as.character(x)
			compl <- (length(which(x %in% c("N", "n", "?", "-", "O", "o", "X", "x"))) / ntax) < propthres
			return(compl)
			})
			
	 if(codon){
			codIDs <- rep(1:(length(compliantSites)/3), each = 3)
			codsToKeep <- rep(F, length(compliantSites))
			for(i in 1:max(codIDs)){
			      if(all(compliantSites[which(codIDs == i)])) codsToKeep[which(codIDs == i)] <- T
			}
			al <- al[, as.logical(codsToKeep)]
			return(al)
	 } else {
			al <- al[, as.logical(compliantSites)]
			return(al)
	 }
}


fasta <- "alignments/Whole_genome_1.fst" ###change this to your input file
aln <- read.dna(fasta, format = "fasta")  
trimmed_seqs <- trimCols(aln, prop = 0.8, codon = F)

#remove gappy chunks at beginning
aln_mat <- as.character(trimmed_seqs)
total_cols <- ncol(aln_mat)
aln_trimmed <- aln_mat[, -(c(1:5081, (total_cols-2060):total_cols))]
output_file <- "alignments/trim_genome_africa.fasta"
write.dna(aln_trimmed, file = output_file, format = "fasta", colsep = "")

#rsync -avz /Users/qjs599/Desktop/mol-evol-phylo/alignments/trim_genome_africa.fasta rachel@ku-01.ifsv-bd.lan.ku.dk:/home/rachel/mol_evol/alignments/trim_genome_africa.fasta