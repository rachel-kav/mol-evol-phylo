
library(Biostrings)
library(ape)
library(seqinr)
library(tidyverse)
library(dplyr)
library(TreeDist)
library(MetBrewer)
library(phytools)
library(cowplot)
setwd("/Users/morah/Library/CloudStorage/OneDrive-ITU/Documents/Github repos/Courses/Molecular Evolution and Phylogenomic course/mol-evol-phylo/mast_africa_output")
trimmed_sitelh <- read.table("mast1.sitelh", header = TRUE)[, 1:5]
colnames(trimmed_sitelh) <- c("site", "LnL", 'aln_segment_132237_to_142408', 'aln_segment_111893_to_122064', 'aln_segment_122065_to_132236')

long_sitelh <- pivot_longer(trimmed_sitelh, cols = c('aln_segment_132237_to_142408', 'aln_segment_111893_to_122064', 'aln_segment_122065_to_132236'), values_to = "likelihood", names_to = "tree")
long_sitelh$LnL <- NULL

long_sitelh <- long_sitelh %>%
  group_by(site) %>%
  mutate(rellh = (max(likelihood)) - (sort(likelihood, decreasing = TRUE)[2]),
         rank = rev(base::rank(likelihood)))

#removing all other ranks than 1, to get just one row / site
rel_lh <- long_sitelh %>% filter(rank==1)

trimmed_sitelh$rel_lh <- rel_lh$rellh


trimmed_sitelh$most_likely <- colnames(trimmed_sitelh[3:5])[apply(trimmed_sitelh[3:5], 1, which.max)]

names(which.max(table(trimmed_sitelh$most_likely)))

window_size = 200

#find the modal class - ie for each 20 bp segment, what is the most likely/probable tree

most_likely_segment <- c()
for (i in 1:dim(trimmed_sitelh)[1]) {
  
  if (i >= dim(trimmed_sitelh)[1]-window_size) {
    segment <- trimmed_sitelh[i:dim(trimmed_sitelh)[1],]
    most_freq <- names(which.max(table(segment$most_likely)))
    most_likely_segment <- append(most_likely_segment, most_freq)
  } else {
    segment <- trimmed_sitelh[i:(i+window_size),]
    most_freq <- names(which.max(table(segment$most_likely)))
    most_likely_segment <- append(most_likely_segment, most_freq)
  }
}

table(most_likely_segment)

trimmed_sitelh$window_lh <- most_likely_segment

modal_class <- c()
modal_rellikelihood <- c()
modal_proportion <- c()

for (i in 1:dim(trimmed_sitelh)[1]) {
  if (i<20) { segment <- trimmed_sitelh[1:(i+(window_size-1)),]
  } else if (i>dim(trimmed_sitelh)[1]-(window_size-1)) {       
    segment <- trimmed_sitelh[(i-(window_size-1)):dim(trimmed_sitelh)[1],]
  } else {
    segment <- trimmed_sitelh[(i-(window_size-1)):(i+(window_size-1)),] 
  }
  most_freq <- names(which.max(table(segment$window_lh)))
  modal_class <- append(modal_class, most_freq)
  modal_rellikelihood <- append(modal_rellikelihood,  mean(segment$rel_lh, na.rm = TRUE))
  modal_proportion <- append(modal_proportion, (sum(segment$window_lh == most_freq)/dim(segment)[1]))
}

table(most_likely_segment)
table(modal_class)
summary(modal_rellikelihood)

trimmed_sitelh$modal_class <- modal_class
trimmed_sitelh$modal_prob <- modal_rellikelihood
trimmed_sitelh$modal_prop <- modal_proportion

trimmed_sitelh$site <- as.numeric(row.names(trimmed_sitelh))

#trimmed_siteprob$modal_class <- factor(trimmed_siteprob$modal_class, levels = c('tree_30', 'tree_8', 'tree_31','tree_10','tree_2',  'tree_21'))

pall <- c("#3c4b99", "#924099", "#df9ed4")




modal_class_plot <- ggplot(trimmed_sitelh, aes(x = site, y = 'tree', colour = modal_class))+
  geom_point(size = 12, shape = 108)+
  theme(axis.title = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        legend.position = 'top')+
  scale_colour_manual(values = pall)+
  scale_x_continuous(expand = c(0, 0), breaks = c(0,10000,20000,30000,40000,50000,60000,70000,80000,90000,100000,110000,120000,130000,140000,150000))#+
#guides(colour = 'none')

modal_class_plot

ggsave("likelihood_modal_class_plot.png", modal_class_plot, device = 'png', width = 380, height = 100, units = 'mm', dpi = 300, bg = 'white')




spline_int <- as.data.frame(spline(trimmed_sitelh$site, trimmed_sitelh$modal_prob))

mean_prob_plot <- ggplot(trimmed_sitelh) + 
  geom_line(data = spline_int, aes(x = x, y = y), alpha = 0.5, color = "grey50")+
  geom_point(aes(x = site, y = modal_prob, colour = modal_class), size = 0.7) +
  scale_x_continuous(expand = c(0, 0), breaks = c(0,10000,20000,30000,40000,50000,60000,70000,80000,90000,100000,110000,120000,130000,140000,150000))+
  scale_y_continuous(breaks = c(0,1,2,3,4))+
  scale_colour_manual(values =pall)+
  theme(axis.title = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = 'bottom')+
  #  labs(y = "Mean probability")+
  guides(colour = 'none')

mean_prob_plot

#save to png 

ggsave("mean_likelihood_modal.png", mean_prob_plot, device = 'png', width = 373, height = 100, units = 'mm', dpi = 400, bg = 'white')




spline_int <- as.data.frame(spline(trimmed_sitelh$site, trimmed_sitelh$modal_prop))

proportion_plot <- ggplot(trimmed_sitelh) + 
  geom_line(data = spline_int, aes(x = x, y = y), color = "grey50", alpha = 0.5)+
  geom_point(aes(x = site, y = modal_prop, colour = modal_class), size = 0.7) +
  scale_x_continuous(expand = c(0, 0), breaks = c(0,10000,20000,30000,40000,50000,60000,70000,80000,90000,100000,110000,120000,130000,140000,150000))+
  scale_colour_manual(values =pall)+
  theme(axis.title = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = 'bottom')+
  labs(y = "Class proportion")+
  guides(colour = 'none')

proportion_plot

#save to png

ggsave("proportion_modal.png", proportion_plot, device = 'png', width = 373, height = 100, units = 'mm', dpi = 400, bg = 'white')

#trimmed_sitelh$windowentropy <- entropydf$rollingmean
#trimmed_sitelh$gappyness <- gapdf$rollingmean

