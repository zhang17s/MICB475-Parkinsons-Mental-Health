#Load libraries
library(tidyverse)
library(phyloseq)
library(DESeq2)
library(cowplot)
library(ggplot2)
library(gridExtra)

#### Load data ####
load("1_make_phyloseq_object/parkinsons_final_depression.RData")
load("1_make_phyloseq_object/parkinsons_final_anxiety.RData")
load("1_make_phyloseq_object/parkinsons_final_sleep.RData")

#### DESeq Volcano Plot for Depression ####
depression_plus1 <- transform_sample_counts(parkinsons_final_depression, function(x) x+1)
depression_deseq <- phyloseq_to_deseq2(depression_plus1, ~`depression_binned`)
DESEQ_depression <- DESeq(depression_deseq)
res_depression <- results(DESEQ_depression, tidy=TRUE, 
                          #this will ensure that No is your reference group
                          contrast = c("depression_binned","Yes","No"))

## Make variable to color by whether it is significant + large change ##
depression_vol_plot <- res_depression %>%
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>%
  ggplot() +
  ggtitle('Depression DeSeq2 Volcano Plot') +
  theme(plot.title = element_text(hjust = 0.5))+
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant))
depression_vol_plot

#### DESeq Volcano Plot for Anxiety ####
anxiety_plus1 <- transform_sample_counts(parkinsons_final_anxiety, function(x) x+1)
anxiety_deseq <- phyloseq_to_deseq2(anxiety_plus1, ~`anxiety_binned`)
DESEQ_anxiety <- DESeq(anxiety_deseq)
res_anxiety <- results(DESEQ_anxiety, tidy=TRUE, 
                       #this will ensure that No is your reference group
                       contrast = c("anxiety_binned","Yes","No"))

## Make variable to color by whether it is significant + large change ##
anxiety_vol_plot <- res_anxiety %>%
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>%
  ggplot() +
  ggtitle('Anxiety DeSeq2 Volcano Plot') +
  theme(plot.title = element_text(hjust = 0.5))+
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant))

#### DESeq Volcano Plot for Sleep ####
sleep_plus1 <- transform_sample_counts(parkinsons_final_sleep, function(x) x+1)
sleep_deseq <- phyloseq_to_deseq2(sleep_plus1, ~`Sleep_problems`)
DESEQ_sleep <- DESeq(sleep_deseq)
res_sleep <- results(DESEQ_sleep, tidy=TRUE, 
                     #this will ensure that No is your reference group
                     contrast = c("Sleep_problems","Yes","No"))

## Make variable to color by whether it is significant + large change ##
sleep_vol_plot <- res_sleep %>%
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>%
  ggplot() +
  ggtitle('Sleep Problem DeSeq2 Volcano Plot') +
  theme(plot.title = element_text(hjust = 0.5))+
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant))

#### Annotate figures with captions “A, B, C, D” ####
combined_volcano <- plot_grid(depression_vol_plot + theme(legend.position="none"), 
                              anxiety_vol_plot + theme(legend.position="none"), 
                              sleep_vol_plot + theme(legend.position="none"),
                              labels = c('A', 'B','C'), nrow = 1)
legend_for_combined_volcano <- get_legend(
  # create some space to the left of the legend
  depression_vol_plot + theme(legend.box.margin = margin(0, 0, 0, 12))
)

combined_volcano_w_legend <-plot_grid(combined_volcano, legend_for_combined_volcano, rel_widths = c(3, .4))
combined_volcano_w_legend

ggsave('4_deseq2/combined_volcano_with_legend.png',
       width = 20,
       height = 10,
       combined_volcano_w_legend)


##### DeSeq2 Bar Plots #####

# To get table of results of depression
sigASVs_depression <- res_depression %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)
View(sigASVs_depression)
# Get only asv names
sigASVs_vec_depression <- sigASVs_depression %>%
  pull(ASV)

# Prune Depression phyloseq file #
depression_DESeq <- prune_taxa(sigASVs_vec_depression,parkinsons_final_depression)
sigASVs_depression <- tax_table(depression_DESeq) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_depression) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))

# generate depression bar plot #
depression_bar_plot <- ggplot(sigASVs_depression) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  ggtitle('Depression DeSeq2 Bar Plot') +
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
depression_bar_plot

# To get table of results of anxiety
sigASVs_anxiety <- res_anxiety %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)
View(sigASVs_anxiety)
# Get only asv names
sigASVs_vec_anxiety <- sigASVs_anxiety %>%
  pull(ASV)

# Prune Anxiety phyloseq file #
anxiety_DESeq <- prune_taxa(sigASVs_vec_anxiety,parkinsons_final_anxiety)
sigASVs_anxiety <- tax_table(anxiety_DESeq) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_anxiety) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))

# generate anxiety bar plot #
anxiety_bar_plot <- ggplot(sigASVs_anxiety) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  ggtitle('Anxiety DeSeq2 Bar Plot') +
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
anxiety_bar_plot

# To get table of results of sleep problems
sigASVs_sleep <- res_sleep %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)
View(sigASVs_sleep)
# Get only asv names
sigASVs_vec_sleep <- sigASVs_sleep %>%
  pull(ASV)

# Prune Sleep Problem phyloseq file #
sleep_problem_DESeq <- prune_taxa(sigASVs_vec_sleep,parkinsons_final_sleep)
sigASVs_sleep <- tax_table(sleep_problem_DESeq) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_sleep) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))

# generate anxiety bar plot #
sleep_problem_bar_plot <- ggplot(sigASVs_sleep) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  ggtitle('Sleep Problem DeSeq2 Bar Plot') +
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
sleep_problem_bar_plot

combined_bar_plot_deseq <- plot_grid(depression_bar_plot, 
                                     anxiety_bar_plot, 
                                     sleep_problem_bar_plot,
                                     labels = c('A', 'B','C'), nrow = 1)
combined_bar_plot_deseq

ggsave('4_deseq2/combined_bar_plot.png',
       width = 20,
       height = 10,
       combined_bar_plot_deseq)
