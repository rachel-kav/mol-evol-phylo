# Requires package TreeDist - Takes the distance between trees where each taxon has been removed as the discordance metric for the taxon.
# The greater the distance, the greater the contribution of a taxon to tree discordance.

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
