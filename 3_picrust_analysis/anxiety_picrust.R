#Install packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

pkgs <- c("phyloseq", "ALDEx2", "SummarizedExperiment", "Biobase", "devtools", 
          "ComplexHeatmap", "BiocGenerics", "BiocManager", "metagenomeSeq", 
          "Maaslin2", "edgeR", "lefser", "limma", "KEGGREST", "DESeq2")

for (pkg in pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE))
    BiocManager::install(pkg)
}

additional_pckgs <- c("ggpicrust2", "tidyverse")
if (any(additional_pckgs== F)) {
  install.packages(packages[!additional_pckgs])
}

#load libraries
library(readr)
library(ggpicrust2)
library(tibble)
library(tidyverse)
library(ggprism)
library(patchwork)
library(DESeq2)
source("3_picrust_analysis/DESeq2_function.R")
library("ggh4x")

#Importing the pathway PICrsut2
abundance_file <- "3_picrust_analysis/Picrust analysis _path_abun_unstrat.tsv"
abundance_data <- read_delim(abundance_file, delim = "\t", col_names = TRUE, trim_ws = TRUE)
abundance_data  =as.data.frame(abundance_data)

#Read in metadata
metadata <- read_delim("1_make_phyloseq_object/parkinsons_metadata_new_edited.csv")

########## PICRUSt on Anxiety Controls ##########
#Filter your metadata as needed to look at specific comparisons
Control_metadata = metadata %>%
  filter(Disease == "Control")

#Remove NAs for anxiety_binned
Control_metadata = Control_metadata[!is.na(Control_metadata$anxiety_binned),]

#Filtering the abundance table to only include samples that are in the filtered metadata
sample_names = Control_metadata$`X.SampleID`
sample_names = append(sample_names, "pathway")
abundance_data_filtered = abundance_data[, colnames(abundance_data) %in% sample_names] #This step is the actual filtering

#Removing individuals with no data that caused a problem for pathways_daa()
abundance_data_filtered =  abundance_data_filtered[, colSums(abundance_data_filtered != 0) > 0]

#Filtering out columns that represent a total abundance < 1000
#Filtering out rows (pathways) that have a total count < 100
abundance_data_filtered = abundance_data_filtered[,colSums(abundance_data_filtered[,-1]) > 1000]
abundance_data_filtered = abundance_data_filtered[rowSums(abundance_data_filtered[,-1]) > 100,]

#Ensuring the rownames for the abundance_data_filtered is empty. This is required for their functions to run.
rownames(abundance_data_filtered) = NULL

#verify samples in metadata match samples in abundance_data
abun_samples = rownames(t(abundance_data_filtered[,-1])) #Getting a list of the sample names in the newly filtered abundance data
Control_metadata = Control_metadata[Control_metadata$`X.SampleID` %in% abun_samples,] #making sure the filtered metadata only includes these samples

#Perform pathway DAA using DESEQ2 method
abundance_daa_results_df <- pathway_daa(abundance = abundance_data_filtered %>% column_to_rownames("pathway"), metadata = Control_metadata, group = "anxiety_binned", daa_method = "DESeq2")

# Annotate MetaCyc pathway results without KO to KEGG conversion
metacyc_daa_annotated_results_df <- pathway_annotation(pathway = "MetaCyc", daa_results_df = abundance_daa_results_df, ko_to_kegg = FALSE)

### Generate pathway heatmap ###
feature_with_p_0.05 <- abundance_daa_results_df %>% filter(p_values < 0.05)

#Changing the pathway column to description for the results 
feature_desc = inner_join(feature_with_p_0.05,metacyc_daa_annotated_results_df, by = "feature")
feature_desc$feature = feature_desc$description
feature_desc = feature_desc[,c(1:7)]
colnames(feature_desc) = colnames(feature_with_p_0.05)

#Changing the pathway column to description for the abundance table
abundance = abundance_data_filtered %>% filter(pathway %in% feature_with_p_0.05$feature)
colnames(abundance)[1] = "feature"
abundance_desc = inner_join(abundance,metacyc_daa_annotated_results_df, by = "feature")
abundance_desc$feature = abundance_desc$description
abundance_desc = abundance_desc[,-c(58:ncol(abundance_desc))]

# Generate pathway PCA plot
anx_control_pca <- pathway_pca(abundance = abundance_data_filtered %>% 
                                 column_to_rownames("pathway"), metadata = Control_metadata, group = "anxiety_binned")
anx_control_pca

#Generate log 2 fold change data for yes vs. no samples
res =  DEseq2_function(abundance_data_filtered,Control_metadata,"anxiety_binned")
res$feature =rownames(res)
res_desc = inner_join(res,metacyc_daa_annotated_results_df, by = "feature")
res_desc = res_desc[, -c(8:13)]
View(res_desc)

#filter to only keep significant pathways
sig_res = res_desc %>%
  filter(pvalue < 0.05)

#make log 2 fold change plot
sig_res <- sig_res[order(sig_res$log2FoldChange),]
anx_control_log <- ggplot(data = sig_res, aes(y = reorder(description, sort(as.numeric(log2FoldChange))), x= log2FoldChange, fill = pvalue))+
  geom_bar(stat = "identity")+ 
  theme_bw()+
  labs(x = "Log2Fold Change", y=" Metabolic Pathway", fill = "P Value") +
  #theme(axis.text.y = element_text(size = 10)) +
  ggtitle("Anxiety Control Cohort") + theme(plot.title=element_text(hjust = 0.5)) +  
  guides(fill = "none") +
  theme(axis.text = element_blank()) +
  theme(axis.title=element_text(size=14,face="bold")) +
  theme(plot.title = element_text(size = 16, face = "bold"))
anx_control_log


#Saving figures for manuscript 
ggsave(filename = "3_picrust_analysis/publication_picrust_figures/fig3_A_control_pca.png", anx_control_pca,
       height = 6, width = 9)

ggsave(filename = "3_picrust_analysis/publication_picrust_figures/fig3_B_control_log.png", anx_control_log, 
       height = 6, width = 10)

########## PICRUSt on Anxiety PD ##########
#Filter for PD patients
PD_metadata = metadata %>%
  filter(Disease == "PD")

#Remove NAs for anxiety
PD_metadata = PD_metadata[!is.na(PD_metadata$anxiety_binned),]

#Filtering the abundance table to only include samples that are in the filtered metadata
sample_names = PD_metadata$`X.SampleID`
sample_names = append(sample_names, "pathway")
abundance_data_filtered = abundance_data[, colnames(abundance_data) %in% sample_names] #This step is the actual filtering

#Removing individuals with no data that caused a problem for pathways_daa()
abundance_data_filtered =  abundance_data_filtered[, colSums(abundance_data_filtered != 0) > 0]

#Filtering out columns that represent a total abundance < 1000
#Filtering out rows (pathways) that have a total count < 100
abundance_data_filtered = abundance_data_filtered[,colSums(abundance_data_filtered[,-1]) > 1000]
abundance_data_filtered = abundance_data_filtered[rowSums(abundance_data_filtered[,-1]) > 100,]

#Ensuring the rownames for the abundance_data_filtered is empty. This is required for their functions to run.
rownames(abundance_data_filtered) = NULL

#verify samples in metadata match samples in abundance_data
abun_samples = rownames(t(abundance_data_filtered[,-1])) #Getting a list of the sample names in the newly filtered abundance data
PD_metadata = PD_metadata[PD_metadata$`X.SampleID` %in% abun_samples,] #making sure the filtered metadata only includes these samples

#Perform pathway DAA using DESeq2 method
abundance_daa_results_df <- pathway_daa(abundance = abundance_data_filtered %>% column_to_rownames("pathway"), metadata = PD_metadata, group = "anxiety_binned", daa_method = "DESeq2")

# Annotate MetaCyc pathway results without KO to KEGG conversion
metacyc_daa_annotated_results_df <- pathway_annotation(pathway = "MetaCyc", daa_results_df = abundance_daa_results_df, ko_to_kegg = FALSE)

### Generate pathway heatmap ###
feature_with_p_0.05 <- abundance_daa_results_df %>% filter(p_values < 0.05)

#Changing the pathway column to description for the results 
feature_desc = inner_join(feature_with_p_0.05,metacyc_daa_annotated_results_df, by = "feature")
feature_desc$feature = feature_desc$description
feature_desc = feature_desc[,c(1:7)]
colnames(feature_desc) = colnames(feature_with_p_0.05)

#Changing the pathway column to description for the abundance table
abundance = abundance_data_filtered %>% filter(pathway %in% feature_with_p_0.05$feature)
colnames(abundance)[1] = "feature"
abundance_desc = inner_join(abundance,metacyc_daa_annotated_results_df, by = "feature")
abundance_desc$feature = abundance_desc$description
abundance_desc = abundance_desc[,-c(77:ncol(abundance_desc))]

# Generate pathway PCA plot
anx_PD_pca <- pathway_pca(abundance = abundance_data_filtered %>% column_to_rownames("pathway"), metadata = PD_metadata, group = "anxiety_binned")
anx_PD_pca

#Generate log 2 fold change for yes. vs no samples
res =  DEseq2_function(abundance_data_filtered,PD_metadata,"anxiety_binned")
res$feature =rownames(res)
res_desc = inner_join(res,metacyc_daa_annotated_results_df, by = "feature")
res_desc = res_desc[, -c(8:13)]
View(res_desc)

#filter to keep only significant pathways
sig_res = res_desc %>%
  filter(pvalue < 0.05)

#Make log 2 fold change plot
sig_res <- sig_res[order(sig_res$log2FoldChange),]
anx_PD_log <- ggplot(data = sig_res, aes(y = reorder(description, sort(as.numeric(log2FoldChange))), x= log2FoldChange, fill = pvalue))+
  geom_bar(stat = "identity")+
  theme_bw()+
  labs(x = "Log2Fold Change", y="Metabolic Pathway", fill = "P Value") +
  #theme(axis.text.y = element_text(size = 10)) +
  ggtitle("Anxiety PD Cohort") + theme(plot.title=element_text(hjust = 0.5)) +
  guides(fill = "none") +
  theme(axis.text = element_blank()) +
  theme(axis.title=element_text(size=14,face="bold")) +
  theme(plot.title = element_text(size = 16, face = "bold"))
anx_PD_log

#saving log plot and pca for manuscript 
ggsave(filename = "3_picrust_analysis/publication_picrust_figures/fig3_C_PD_pca.png", anx_PD_pca,
       height = 6, width = 9)

ggsave(filename = "3_picrust_analysis/publication_picrust_figures/fig3_D_PD_log.png", anx_PD_log, 
       height = 6, width = 15)



