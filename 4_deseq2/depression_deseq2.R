#Load libraries
library(tidyverse)
library(phyloseq)
library(DESeq2)

#### Load data ####
load("1_make_phyloseq_object/parkinsons_final_depression.RData")

#### DESeq for Depression ####
depression_plus1 <- transform_sample_counts(parkinsons_final_depression, function(x) x+1)
depression_deseq <- phyloseq_to_deseq2(depression_plus1, ~`depression_binned`)
DESEQ_depression <- DESeq(depression_deseq)
res <- results(DESEQ_depression, tidy=TRUE, 
               #this will ensure that No is your reference group
               contrast = c("depression_binned","Yes","No"))
View(res)

## Make variable to color by whether it is significant + large change ##
depression_vol_plot <- res %>%
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant)) +
  ggtitle("Depression DeSeq2 Volcano Plot") + theme(plot.title=element_text(hjust=0.5))

ggsave(filename="depression_volcano_plot.png",depression_vol_plot)

# To get table of results
sigASVs <- res %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)
View(sigASVs)

# Get only asv names
sigASVs_vec <- sigASVs %>%
  pull(ASV)

# Prune phyloseq file
depression_DESeq <- prune_taxa(sigASVs_vec,parkinsons_final_depression)
sigASVs <- tax_table(depression_DESeq) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))

depression_bar_plot <- ggplot(sigASVs) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
depression_bar_plot

ggsave('4_deseq2/depression_bar_plot.png',
       depression_bar_plot)
