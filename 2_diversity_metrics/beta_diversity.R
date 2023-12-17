#Load libraries
library(phyloseq)
library(ape) 
library(tidyverse)
library(vegan)
library(cowplot)
library(gridExtra)

### load rarefied phyloseq object ###
load("1_make_phyloseq_object/parkinsons_edited_rare.RData")

#Filter out PD patients 
PD_patients <- subset_samples(parkinsons_final_anxiety, `Disease` == "PD")

#Filter out control
Ctrl_patients <- subset_samples(parkinsons_final_anxiety, `Disease` == "Control")

###Data frame 
#PD
samp_dat_wdiv_PD <- data.frame(sample_data(PD_patients), estimate_richness(PD_patients))

#Control
samp_dat_wdiv_Ctrl <- data.frame(sample_data(Ctrl_patients), estimate_richness(Ctrl_patients))

## Jaccard ## 
#PD patients
jac_dm <- distance(parkinsons_rare, method = "jaccard", binary = T)
pcoa_jac_PD <- ordinate(parkinsons_rare, method = "NMDS", distance = jac_dm)
PD_jac <- plot_ordination(parkinsons_rare, pcoa_jac_PD, color = "Disease") +
  labs(col = "Disease Status") + theme_bw() + stat_ellipse(level = 0.95) +
  ggtitle("Jaccard Disease") + theme(plot.title = element_text(hjust = 0.5)) +
  guides(fill = "none") 
PD_jac <- plot_grid(PD_jac, labels = c('G'))

ggsave("PD_jac_pcoa.png"
       , PD_jac
       , height=4, width=6)
PD_jac



### Beta Diversity ###
## Jaccard ## 
jac_dm <- distance(parkinsons_rare, method = "jaccard", binary = TRUE)
pcoa_jac <- ordinate(parkinsons_rare, method = "NMDS", distance = jac_dm)
gg_jac_pcoa <- plot_ordination(parkinsons_rare, pcoa_jac, color = "Disease") +
  labs(col = "Disease Status") + theme_bw() + stat_ellipse(level = 0.95) +
  ggtitle("Jaccard_Disease") + theme(plot.title = element_text(hjust = 0.5)) +
  ggside::geom_xsideboxplot(aes(fill = Disease, y = Disease), 
                            orientation = "y") +
  ggside::geom_ysideboxplot(aes(fill = Disease, x = Disease), 
                            orientation = "x") +
  ggside::scale_xsidey_discrete(labels = NULL) +
  ggside::scale_ysidex_discrete(labels = NULL) +
  ggside::theme_ggside_void()

gg_jac_pcoa

ggsave("Disease_jac_pcoa.png"
       , gg_jac_pcoa
       , height=4, width=6)

adonis2(jac_dm ~ Disease, data = samp_dat_wdiv)

## bray curtis ##
bray_dm <- distance(parkinsons_rare, method = "bray")
pcoa_bray <- ordinate(parkinsons_rare, method = "PCoA", distance = bray_dm)
gg_bray_pcoa <- plot_ordination(parkinsons_rare, pcoa_bray, color = "Disease") +
  labs(col = "Disease Status") + theme_bw() + stat_ellipse(level = 0.95) +
  ggtitle("Bray_Curtis_Disease") + theme(plot.title = element_text(hjust = 0.5)) +
  ggside::geom_xsideboxplot(aes(fill = Disease, y = Disease), 
                            orientation = "y") +
  ggside::geom_ysideboxplot(aes(fill = Disease, x = Disease), 
                            orientation = "x") +
  ggside::scale_xsidey_discrete(labels = NULL) +
  ggside::scale_ysidex_discrete(labels = NULL) +
  ggside::theme_ggside_void() 
gg_bray_pcoa

ggsave("Disease_bray_pcoa.png"
       , gg_bray_pcoa
       , height=4, width=6)

adonis2(bray_dm ~ Disease, data = samp_dat_wdiv)

## unweighted unifrac ##
unifrac_dm <- distance(parkinsons_rare, method = "unifrac")
pcoa_unifrac <- ordinate(parkinsons_rare, method = "PCoA", distance = unifrac_dm)
gg_unifrac_pcoa <- plot_ordination(parkinsons_rare, pcoa_unifrac, color = "Disease") +
  labs(col = "Disease Status") + theme_bw() + stat_ellipse(level = 0.95) +
  ggtitle("Unweighted_Unifrac_Disease") + theme(plot.title = element_text(hjust = 0.5))+
  ggside::geom_xsideboxplot(aes(fill = Disease, y = Disease), 
                            orientation = "y") +
  ggside::geom_ysideboxplot(aes(fill = Disease, x = Disease), 
                            orientation = "x") +
  ggside::scale_xsidey_discrete(labels = NULL) +
  ggside::scale_ysidex_discrete(labels = NULL) +
  ggside::theme_ggside_void() 
gg_unifrac_pcoa

ggsave("Disease_unifrac_pcoa.png"
       , gg_unifrac_pcoa
       , height=4, width=6)

adonis2(unifrac_dm ~ Disease, data = samp_dat_wdiv)

## weighted_unifrac ##
w_unifrac_dm <- distance(parkinsons_rare, method ="wunifrac")
pcoa_w_unifrac <- ordinate(parkinsons_rare, method="PCoA", distance=w_unifrac_dm)
gg_wunifrac_pcoa <- plot_ordination(parkinsons_rare, pcoa_w_unifrac, color = "Disease") +
  labs(col = "Disease Status") + theme_bw() + stat_ellipse(level = 0.95) +
  ggtitle("Weighted_Unifrac_Disease") + theme(plot.title = element_text(hjust = 0.5)) +
  ggside::geom_xsideboxplot(aes(fill = Disease, y = Disease), 
                            orientation = "y") +
  ggside::geom_ysideboxplot(aes(fill = Disease, x = Disease), 
                            orientation = "x") +
  ggside::scale_xsidey_discrete(labels = NULL) +
  ggside::scale_ysidex_discrete(labels = NULL) +
  ggside::theme_ggside_void() 
gg_wunifrac_pcoa

ggsave("Disease_wunifrac_pcoa.png"
       , gg_wunifrac_pcoa
       , height=4, width=6)

adonis2(w_unifrac_dm ~ Disease, data = samp_dat_wdiv)