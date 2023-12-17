#Load libraries
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

gg_richness_PD_anxiety <- plot_richness(PD_patients, x = "anxiety_binned", measures = c("shannon")) + 
  xlab("PD Anxiety") + geom_boxplot()
gg_richness_Ctrl_anxiety <- plot_richness(Ctrl_patients, x = "anxiety_binned", measures = c("shannon")) + 
  xlab("Control Anxiety") + geom_boxplot()

anxiety_shannon <- plot_grid(gg_richness_PD_anxiety, gg_richness_Ctrl_anxiety, labels = c('C', 'D'))
anxiety_shannon

### Load depression phyloseq Object 
load("1_make_phyloseq_object/parkinsons_final_depression.RData")

#Filter out PD patients 
PD_patients <- subset_samples(parkinsons_final_depression, `Disease` == "PD")

#Filter out control
Ctrl_patients <- subset_samples(parkinsons_final_depression, `Disease` == "Control")

gg_richness_PD_depression <- plot_richness(PD_patients, x = "depression_binned", measures = c("shannon")) + 
  xlab("PD Depression") + geom_boxplot()
gg_richness_Ctrl_depression <- plot_richness(Ctrl_patients, x = "depression_binned", measures = c("shannon")) + 
  xlab("Control Depression") + geom_boxplot()

depression_shannon <- plot_grid(gg_richness_PD_depression, gg_richness_Ctrl_depression, labels = c('A', 'B'))
depression_shannon

### Load depression phyloseq Object 
load("1_make_phyloseq_object/parkinsons_final_sleep_problems.RData")

#Filter out PD patients 
PD_patients <- subset_samples(parkinsons_final_sleep_problems, `Disease` == "PD")

#Filter out control
Ctrl_patients <- subset_samples(parkinsons_final_sleep_problems, `Disease` == "Control")

gg_richness_PD_sleep <- plot_richness(PD_patients, x = "Sleep_problems", measures = c("shannon")) + 
  xlab("PD Sleep Problems") + geom_boxplot()
gg_richness_Ctrl_sleep <- plot_richness(Ctrl_patients, x = "Sleep_problems", measures = c("shannon")) + 
  xlab("Control Sleep Problems") + geom_boxplot()

sleep_shannon <- plot_grid(gg_richness_PD_sleep, gg_richness_Ctrl_sleep, labels = c('E', 'F'))
sleep_shannon

### load rarefied phyloseq object ###
load("1_make_phyloseq_object/parkinsons_edited_rare.RData")

gg_richness_PD <- plot_richness(parkinsons_rare, x = "Disease", measures = c("shannon")) + 
  xlab("Disease Status") + geom_boxplot()

PD_shannon <- plot_grid(gg_richness_PD,  labels = c('G'))
PD_shannon

dep_anxiety_sleep_disease_together_shannon <- grid.arrange(depression_shannon, anxiety_shannon,sleep_shannon, PD_shannon, ncol = 1)
dep_anxiety_sleep_disease_together_shannon

ggsave("2_diversity_metrics/Supplementary_Fig2_Shannon.png"
       , dep_anxiety_sleep_disease_together_shannon
       , height=15, width=10)
