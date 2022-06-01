# Script for preparing the RNAseq data to upload to CIBERSORTx

#Setup environment
rm(list = ls())

# Install packages and load
required.bioconductor.packages = c("GSVA","ComplexHeatmap", "gclus")                                                                   
library(required.bioconductor.packages)

# Set parameters
#cluster = "S6"
#Cancer = "CCSK"
# Load data
load("./Processed_Data/003_PanCancer_pediatric_only_normalized_TCGAbiolinks_TP_filtered.Rdata")
#load("./Analysis/ICR_data_TARGET_PanCancer/ICR_clustered_PanCancers_after_split/TARGET_AllCancers_together_table_cluster_with_Cancername.RData")
load("./Analysis/after_split/Signature_Enrichment/Thorsson_clustering/030.Thorsson_samples_cancers_10000repeats_HM_clusters_km6_pediatric_only_V3.Rdata")
#ICR_all = ICR_all[which(ICR_all$Cancer %in% Cancer),]
#rownames(ICR_all) = gsub("\\.","-",rownames(ICR_all))
#annotation = annotation[which(annotation$cluster %in% cluster),]

#filtered.norm.RNAseqData = filtered.norm.RNAseqData[,which(colnames(filtered.norm.RNAseqData) %in% annotation$Sample)] 
#load(paste0("./Processed_Data/TARGET_Pan_cancer/",cluster,"_PanNormalized_matrix.Rdata"))
# Instructions for preparation CIBERSORTx
# Tab-delimited tabular input format (.txt) with no double quotations and no missing entries.
# Mixture file formatting requirements. Improperly formatted files may cause CIBERSORTx to fail to run.

#1. Gene symbols in column 1; Mixture labels in row 1
#2. Given the significant difference between counts (e.g., CPM) and gene length-normalized expression data (e.g., TPM) we recommend that the signature matrix and mixture files be represented in the same normalization space whenever possible.
#3. Data should be in non-log space. Note: if maximum expression value is less than 50; CIBERSORTx will assume that data are in log space, and will anti-log all expression values by 2x.
#4. CIBERSORTx will add an unique identifier to each redundant gene symbol, however we recommend that users remove redundancy prior to file upload.
#5. CIBERSORTx performs a feature selection and therefore typically does not use all genes in the signature matrix. It is generally ok if some genes are missing from the user’s mixture file. If less than 50% of signature matrix genes overlap, CIBERSORTx will issue a warning.

Mixture_file = filtered.norm.RNAseqData # Data should not be in log-scale
max(Mixture_file) # Maximum value is > 50, so CIBERSORTx will not assume data is in log-scale
rownames(Mixture_file) # rownames are genes (rownames are unique, so no redundant gene symbols (see point 4))
colnames(Mixture_file) # colnames are samples

Mixture_file = data.frame(Mixture_file)
Mixture_file$GeneSymbol = rownames(Mixture_file)
Mixture_file = Mixture_file[, c(ncol(Mixture_file), 1: ncol(Mixture_file) - 1)]

rownames(Mixture_file) = NULL
#colnames(Mixture_file) = gsub("X", "", colnames(Mixture_file))

dir.create(paste0("./Analysis/after_split"), showWarnings = FALSE)
dir.create(paste0("./Analysis/after_split/CIBERSORTx"), showWarnings = FALSE)
#dir.create(paste0("./Analysis/after_split/CIBERSORTx/",cluster,"/"), showWarnings = FALSE)


write.table(Mixture_file, file = paste0("./Analysis/after_split/CIBERSORTx/052.all_408_patients.txt"), sep = "\t")

#######################################################################################
#After analysis by CYBERSORTx

#Calculate and combine the stromal scores for all signatures 
# Setup environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(master.location)
setwd(paste0(master.location,"/TBI-LAB - Project - Pediatric Pan Cancer TARGET"))

S1 = read.csv(paste0("./Analysis/after_split/CIBERSORTx/S1/CIBERSORTx__Thorsson_S1.csv"))

S2 = read.csv(paste0("./Analysis/after_split/CIBERSORTx/S2/CIBERSORTx__Thorsson_S2.csv"))

S3 = read.csv(paste0("./Analysis/after_split/CIBERSORTx/S3/CIBERSORTx__Thorsson_S3.csv"))

S4 = read.csv(paste0("./Analysis/after_split/CIBERSORTx/S4/CIBERSORTx__Thorsson_S4.csv"))

S5 = read.csv(paste0("./Analysis/after_split/CIBERSORTx/S5/CIBERSORTx__Thorsson_S5.csv"))

S6 = read.csv(paste0("./Analysis/after_split/CIBERSORTx/S6/CIBERSORTx__Thorsson_S6.csv"))

Cyberall_all = rbind(S1,S2,S3,S4,S5,S6)

Cyberall_all$Mixture = gsub("\\.","-",Cyberall_all$Mixture)

save(Cyberall_all,file = "./Analysis/after_split/CIBERSORTx/Cybersortx_all_clusters_pediatric_only.Rdata")
