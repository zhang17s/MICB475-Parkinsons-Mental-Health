############ Anxiety ###########
### Load anxiety phyloseq Object 
load("1_make_phyloseq_objects/parkinsons_final_anxiety.RData")

#Filter out PD patients 
PD_patients <- subset_samples(parkinsons_final_anxiety, `Disease` == "PD")

#Filter out control
Ctrl_patients <- subset_samples(parkinsons_final_anxiety, `Disease` == "Control")

#Make data frame for PD
samp_dat_wdiv_PD <- data.frame(sample_data(PD_patients), estimate_richness(PD_patients))

#Make data frame for control
samp_dat_wdiv_Ctrl <- data.frame(sample_data(Ctrl_patients), estimate_richness(Ctrl_patients))

## Jaccard ## 
#PD patients
jac_dm <- distance(PD_patients, method = "jaccard", binary = T)
pcoa_jac_PD <- ordinate(PD_patients, method = "NMDS", distance = jac_dm)
PD_anxiety_jac <- plot_ordination(PD_patients, pcoa_jac_PD, color = "anxiety_binned") +
  labs(col = "") + theme_bw() + stat_ellipse(level = 0.95) +
  ggtitle("Jaccard PD Anxiety") + theme(plot.title = element_text(hjust = 0.5)) 
PD_anxiety_jac

ggsave("publication_figures/diversity/PD_anxiety_jac_pcoa.png"
       , PD_anxiety_jac
       , height=4, width=6)

#Statistical analysis for table 1 data
adonis2(jac_dm ~ `anxiety_binned`, data = samp_dat_wdiv_PD)

#Healthy Controls
jac_dm_ctrl <- distance(Ctrl_patients, method = "jaccard", binary = T)
pcoa_jac_ctrl <- ordinate(Ctrl_patients, method = "NMDS", distance = jac_dm_ctrl)
ctrl_anxiety_jac <- plot_ordination(Ctrl_patients, pcoa_jac_ctrl, color = "anxiety_binned") + 
  labs(col = "Anxiety") + theme_bw() + stat_ellipse(level = 0.95) +
  ggtitle("Jaccard Control Anxiety") + theme(plot.title = element_text(hjust = 0.5)) 
ctrl_anxiety_jac

ggsave("publication_figures/diversity/ctrl_anxiety_jac_pcoa.png"
       , ctrl_anxiety_jac
       , height=4, width=6)

#statistical analysis for table 1 data
adonis2(jac_dm_ctrl ~ `anxiety_binned`, data = samp_dat_wdiv_Ctrl)

anxiety_jac <- plot_grid(PD_anxiety_jac, ctrl_anxiety_jac, labels = c('C', 'D'))
anxiety_jac

################ Depression ################
### Load depression phyloseq Object 
load("1_make_phyloseq_objects/parkinsons_final_depression.RData")

#Filter out PD patients 
PD_patients <- subset_samples(parkinsons_final_depression, `Disease` == "PD")

#Filter out control
Ctrl_patients <- subset_samples(parkinsons_final_depression, `Disease` == "Control")

# Make data frame for PD
samp_dat_wdiv_PD <- data.frame(sample_data(PD_patients), estimate_richness(PD_patients))

#make data frame for Control
samp_dat_wdiv_Ctrl <- data.frame(sample_data(Ctrl_patients), estimate_richness(Ctrl_patients))

## Jaccard ## 
#PD patients
jac_dm <- distance(PD_patients, method = "jaccard", binary = T)
pcoa_jac_PD <- ordinate(PD_patients, method = "NMDS", distance = jac_dm)
PD_depression_jac <- plot_ordination(PD_patients, pcoa_jac_PD, color = "depression_binned") +
  labs(col = "") + theme_bw() + stat_ellipse(level = 0.95) +
  ggtitle("Jaccard PD Depression") + theme(plot.title = element_text(hjust = 0.5)) +
  guides(fill = "none") 
PD_depression_jac

ggsave("publication_figures/diversity/PD_depression_jac_pcoa.png"
       , PD_depression_jac
       , height=4, width=6)

#statistical analysis for table 1 data
adonis2(jac_dm ~ `depression_binned`, data = samp_dat_wdiv_PD)

#Healthy Controls
jac_dm_ctrl <- distance(Ctrl_patients, method = "jaccard", binary = T)
pcoa_jac_ctrl <- ordinate(Ctrl_patients, method = "NMDS", distance = jac_dm_ctrl)
ctrl_depression_jac <- plot_ordination(Ctrl_patients, pcoa_jac_ctrl, color = "depression_binned") + 
  labs(col = "Depression") + theme_bw() + stat_ellipse(level = 0.95) +
  ggtitle("Jaccard Control Depression") + theme(plot.title = element_text(hjust = 0.5)) 
ctrl_depression_jac

ggsave("publication_figures/diversity/ctrl_depression_jac_pcoa.png"
       , ctrl_depression_jac
       , height=4, width=6)

#statistical analysis for table 1 data
adonis2(jac_dm_ctrl ~ `depression_binned`, data = samp_dat_wdiv_Ctrl)

depression_jac <- plot_grid(PD_depression_jac, ctrl_depression_jac, labels = c('A', 'B'))
depression_jac

dep_anxiety_together <- grid.arrange(depression_jac, anxiety_jac, ncol = 1)
dep_anxiety_together