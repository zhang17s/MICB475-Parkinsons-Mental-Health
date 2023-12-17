#load libraries
library(phyloseq)
library(ape) 
library(tidyverse)
library(vegan)
library(cowplot)
library(gridExtra)

############ Anxiety ###########
### Load anxiety phyloseq Object 
load("1_make_phyloseq_object/parkinsons_final_anxiety.RData")

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
jac_dm <- distance(PD_patients, method = "jaccard", binary = T)
pcoa_jac_PD <- ordinate(PD_patients, method = "NMDS", distance = jac_dm)
PD_anxiety_jac <- plot_ordination(PD_patients, pcoa_jac_PD, color = "anxiety_binned") +
  labs(col = "") + theme_bw() + stat_ellipse(level = 0.95) +
  ggtitle("Jaccard PD Anxiety") + theme(plot.title = element_text(hjust = 0.5)) 
PD_anxiety_jac

#Stats analysis for table 1
adonis2(jac_dm ~ `anxiety_binned`, data = samp_dat_wdiv_PD)

#Healthy Controls
jac_dm_ctrl <- distance(Ctrl_patients, method = "jaccard", binary = T)
pcoa_jac_ctrl <- ordinate(Ctrl_patients, method = "NMDS", distance = jac_dm_ctrl)
ctrl_anxiety_jac <- plot_ordination(Ctrl_patients, pcoa_jac_ctrl, color = "anxiety_binned") + 
  labs(col = "Anxiety") + theme_bw() + stat_ellipse(level = 0.95) +
  ggtitle("Jaccard Control Anxiety") + theme(plot.title = element_text(hjust = 0.5)) 
ctrl_anxiety_jac

#Stats analysis for table 1
adonis2(jac_dm_ctrl ~ `anxiety_binned`, data = samp_dat_wdiv_Ctrl)

#Label PD and control anxiety jaccard metrics
anxiety_jac <- plot_grid(PD_anxiety_jac, ctrl_anxiety_jac, labels = c('C', 'D'))
anxiety_jac

################ Depression ################
### Load depression phyloseq Object 
load("1_make_phyloseq_object/parkinsons_final_depression.RData")

#Filter out PD patients 
PD_patients <- subset_samples(parkinsons_final_depression, `Disease` == "PD")

#Filter out control
Ctrl_patients <- subset_samples(parkinsons_final_depression, `Disease` == "Control")

###Data frame 
#PD
samp_dat_wdiv_PD <- data.frame(sample_data(PD_patients), estimate_richness(PD_patients))

#Control
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

#Stats analysis for table 1
adonis2(jac_dm ~ `depression_binned`, data = samp_dat_wdiv_PD)

#Healthy Controls
jac_dm_ctrl <- distance(Ctrl_patients, method = "jaccard", binary = T)
pcoa_jac_ctrl <- ordinate(Ctrl_patients, method = "NMDS", distance = jac_dm_ctrl)
ctrl_depression_jac <- plot_ordination(Ctrl_patients, pcoa_jac_ctrl, color = "depression_binned") + 
  labs(col = "Depression") + theme_bw() + stat_ellipse(level = 0.95) +
  ggtitle("Jaccard Control Depression") + theme(plot.title = element_text(hjust = 0.5)) 
ctrl_depression_jac

#Stats analysis for table 1
adonis2(jac_dm_ctrl ~ `depression_binned`, data = samp_dat_wdiv_Ctrl)

#Label PD and control depression metrics 
depression_jac <- plot_grid(PD_depression_jac, ctrl_depression_jac, labels = c('A', 'B'))
depression_jac

######## Sleep problems ########
load("1_make_phyloseq_object/parkinsons_final_sleep.RData")
parkinsons_final_sleep <- subset_samples(parkinsons_final_sleep, !is.na(Sleep_problems))

#Filter out PD patients 
PD_patients <- subset_samples(parkinsons_final_sleep_problems , `Disease` == "PD")

#Filter out control
Ctrl_patients <- subset_samples(parkinsons_final_sleep_problems , `Disease` == "Control")

###Data frame 
#PD
samp_dat_wdiv_PD <- data.frame(sample_data(PD_patients), estimate_richness(PD_patients))

#Control
samp_dat_wdiv_Ctrl <- data.frame(sample_data(Ctrl_patients), estimate_richness(Ctrl_patients))

## Jaccard ## 
#PD patients
jac_dm <- distance(PD_patients, method = "jaccard", binary = T)
pcoa_jac_PD <- ordinate(PD_patients, method = "NMDS", distance = jac_dm)
PD_Sleep_problems_jac <- plot_ordination(PD_patients, pcoa_jac_PD, color = "Sleep_problems") +
  labs(col = "") + theme_bw() + stat_ellipse(level = 0.95) +
  ggtitle("Jaccard PD Sleep") + theme(plot.title = element_text(hjust = 0.5)) +
  guides(fill = "none") 
PD_Sleep_problems_jac

#Stats analysis for table 1
adonis2(jac_dm ~ `Sleep_problems`, data = samp_dat_wdiv_PD)

#Healthy Controls
jac_dm_ctrl <- distance(Ctrl_patients, method = "jaccard", binary = T)
pcoa_jac_ctrl <- ordinate(Ctrl_patients, method = "NMDS", distance = jac_dm_ctrl)
ctrl_Sleep_problems_jac <- plot_ordination(Ctrl_patients, pcoa_jac_ctrl, color = "Sleep_problems") + 
  labs(col = "Sleep problems") + theme_bw() + stat_ellipse(level = 0.95) +
  ggtitle("Jaccard Control Sleep") + theme(plot.title = element_text(hjust = 0.5)) 
ctrl_Sleep_problems_jac

#Stats analysis for table 1
adonis2(jac_dm_ctrl ~ `Sleep_problems`, data = samp_dat_wdiv_Ctrl)

Sleep_problems_jac <- plot_grid(PD_Sleep_problems_jac, ctrl_Sleep_problems_jac, labels = c('E', 'F'))
Sleep_problems_jac
##### PD #####

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

dep_anxiety_sleep_disease_together <- grid.arrange(depression_jac, anxiety_jac,Sleep_problems_jac, PD_jac, ncol = 1)
dep_anxiety_sleep_disease_together
ggsave("Jac_Pcoa.png"
       , dep_anxiety_sleep_disease_together
       , height=15, width=10)