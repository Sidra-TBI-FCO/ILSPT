# Script for the Deconvolution and GSEA

# Setup environment
rm(list=ls())

required.bioconductor.packages = c("GSVA","gclus")                                                                   
library(required.bioconductor.packages)

# Set parameter
Cancer = "PanCancer"
Gene.set = "Bindea_REV1"  # alternative options: "Bindea_REV1" or "selected.pathways"
# Load data and R scripts
load(paste0("~/bioinformatics tools/GSEA list/immune.gene.lists.v4.Rdata"))
load(paste0("./Processed_Data/003_",Cancer,"_pediatric_only_normalized_TCGAbiolinks_TP_filtered.Rdata"))

# Create folders 
dir.create("./Figures/after_split/",showWarnings = FALSE)
dir.create(paste0("./Figures/after_split/Heatmaps/", Gene.set), showWarnings = FALSE)

Expression.data = log(filtered.norm.RNAseqData +1, 2)
available_genes = rownames(Expression.data)
Gene.list = Angelove_ORIG
unavailable_genes_RNAseq = unlist(Gene.list)[-which(unlist(Gene.list) %in% rownames(Expression.data))]
  
cat(paste0(Gene.set," ssGSEA ", ". Total number of genes is ", length(unlist(Gene.list)), ".",
             " Of which ", length(unlist(Gene.list)[unlist(Gene.list) %in% available_genes]), 
             " genes are available in expression data."), append = TRUE, sep = "\n")
  
## ssGSEA
  ES = gsva(Expression.data,Gene.list,method="ssgsea")
 ESz = ES 
  for(j in 1: nrow(ESz))  {
  ESz[j,] = (ES[j,]-mean(ES[j,]))/sd(ES[j,]) # z-score the enrichment matrix
  }
 
  # create directories
  dir.create("./Analysis/after_split/",showWarnings = FALSE)
  dir.create("./Analysis/after_split/Signature_Enrichment", showWarnings = FALSE)
  
  ## Save Scores
  save(ES,ESz, 
       file = paste0("./Analysis/after_split/Signature_Enrichment/GSEA_", Cancer, 
                     "_pediatric_only_", Gene.set,".Rdata"))
