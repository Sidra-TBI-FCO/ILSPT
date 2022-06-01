# Setup environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(master.location)
setwd(paste0(master.location,"/TBI-LAB - Project - Pediatric Pan Cancer TARGET"))
source(paste0(toolbox.path,"/R scripts/ipak.function.R"))
required.packages <- c("corrplot", "stringr")
ipak(required.packages)    

# Set Parameters
Cancer = "OS"  #NBL_mycn_Amp  #NBL_mycn_NA #NBL_mycn_Namp_Intermed.Low  #NBL_mycn_Namp_High
Cancer_skip = c("")
colpattern = colorRampPalette(c("#152B7E", "white", "#1B7E09"))(n = 410)
selected_genes = "Selected.pathways"                                                                                                   # Specify which genes will be correlated
test = "pearson"
score = "ICR_score"

load("./Analysis/after_split/Signature_Enrichment/Thorsson_clustering/030.Thorsson_samples_cancers_10000repeats_HM_clusters_km6_pediatric_only_V3.Rdata")
load(paste0("./Analysis/after_split/Signature_Enrichment/GSEA_PanCancer_pediatric_onlySelected.pathways.Rdata"))
annotation =annotation_all
annotation = annotation[which(annotation$Cancer %in% Cancer),]
ES = ES[,which(colnames(ES) %in% annotation$Sample)]

# Create folders
dir.create("./Figures/after_split/",showWarnings = FALSE)
dir.create(paste0("./Figures/after_split/Correlation_plots/",score), showWarnings = FALSE)
dir.create(paste0("./Figures/after_split/Correlation_plots/",score,"_",selected_genes, "_",score,"_Correlation_plots"), showWarnings = FALSE)

Hallmark.enrichment.score.df = as.data.frame(t(ES))

####### Calculating the scores #######
load("./Processed_Data/003_PanCancer_pediatric_only_normalized_TCGAbiolinks_TP_filtered.Rdata")

filtered.norm.RNAseqData =filtered.norm.RNAseqData[,which(colnames(filtered.norm.RNAseqData) %in% annotation$Sample)]

if (score == "Cytolytic_score") {
  signature1 = c("GZMA", "PRF1", "GNLY", "GZMH", "GZMM")
}

if (score == "ICR_act_score") {
  signature1 = c("IFNG","IRF1", "STAT1","IL12B","TBX21","GZMB","GNLY","PRF1","GZMH", "GZMA","CXCL10","CXCL9","CCL5","CD8A","CD8B","PDCD1")
}

if (score == "ICR_inh_score") {
  signature1 = c("CTLA4", "FOXP3","IDO1" , "CD274")
}

if (score == "ICR_score") {
  signature1 = c("IFNG","IRF1", "STAT1","IL12B","TBX21","GZMB","GNLY","PRF1","GZMH", "GZMA","CXCL10","CXCL9","CCL5","CD8A","CD8B","PDCD1","CTLA4", "FOXP3","IDO1" , "CD274")
}


filtered.norm.RNAseqData = filtered.norm.RNAseqData[which(rownames(filtered.norm.RNAseqData) %in% signature1),]

subset_RNASeq = t(filtered.norm.RNAseqData)
selected_subset_RNAseq_log2 = log(subset_RNASeq +1, 2)
selected_subset_RNAseq_log2 = as.data.frame(selected_subset_RNAseq_log2)
selected_subset_RNAseq_log2[,score] = rowMeans(selected_subset_RNAseq_log2) 
###

Hallmark.enrichment.score.df[,score] = selected_subset_RNAseq_log2[,score][match(row.names(Hallmark.enrichment.score.df), row.names(selected_subset_RNAseq_log2))]

Hallmark_GSEA_cor <- cor (Hallmark.enrichment.score.df,method=test)

mean_correlation_table = data.frame(Cancers = Cancer, Mean.correlation = 0)

cor.mtest <- function(mat, conf.level = 0.95) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat <- lowCI.mat <- uppCI.mat <- matrix(NA, n, n)
  diag(p.mat) <- 0
  diag(lowCI.mat) <- diag(uppCI.mat) <- 1
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], conf.level = conf.level)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
      lowCI.mat[i, j] <- lowCI.mat[j, i] <- tmp$conf.int[1]
      uppCI.mat[i, j] <- uppCI.mat[j, i] <- tmp$conf.int[2]
    }
  }
  return(list(p.mat, lowCI.mat, uppCI.mat))
}
Hallmark_GSEA_cor_sign <- cor.mtest(Hallmark_GSEA_cor, 0.95)

dir.create(paste0("./Analysis/after_split/Correlations/",selected_genes),showWarnings = FALSE)
dir.create(paste0("./Analysis/after_split/Correlations/",selected_genes,"/",score),showWarnings = FALSE)

save(Hallmark_GSEA_cor, Hallmark_GSEA_cor_sign, file = paste0("./Analysis/after_split/Correlations/",selected_genes,"/",score,"/",selected_genes, "_",Cancer,"_Correlation_plots.Rdata"))

###########################################

rm(list=ls())
load("~/R.Config.Rdata")
setwd(master.location)
setwd(paste0(master.location,"/TBI-LAB - Project - Pediatric Pan Cancer TARGET"))
# Install packages and load
source(paste0(toolbox.path,"/R scripts/ipak.function.R"))

required.packages <- c("corrplot", "stringr","ComplexHeatmap","ggplot2")
ipak(required.packages)                                                                                                  # Install and load required packages

# Set Parameters
CancerTYPES = "ALL"                                                                                                    # Specify download method (this information to be used when saving the file)
Cancer_skip = c("")
display_correlations = "irrespective_of_significance"                                                                    # Can either be "only_significant" or "irrespective_of_significance"
my.palette = colorRampPalette(c("#152B7E", "white", "#1B7E09"))(n = 444)

selected_genes = "Selected.Pathways"                                                                                                   # Specify which genes will be correlated
test = "pearson"
Cancer = "Pancancer"
Gene.Set = "Selected.Pathways"
score = "ICR_score"
# Load data
TARGET.cancersets = read.csv(paste0("./Data/TARGET.dataset.csv"),stringsAsFactors = FALSE)                    # TARGET.datasets.csv  (Cancer Types Abbreviations) 
TARGET.cancersets = TARGET.cancersets[-c(1,3,4,5,7),] 

if (CancerTYPES == "ALL") { 
  CancerTYPES = c("Pancancer", TARGET.cancersets$cancerType)
}

N.sets = length(CancerTYPES)

load("./Analysis/after_split/Signature_Enrichment/Thorsson_clustering/030.Thorsson_samples_cancers_10000repeats_HM_clusters_km6_pediatric_only_V3.Rdata")
load(paste0("./Analysis/after_split/Correlations/",selected_genes,"/",score,"/Selected.pathways_CCSK_Correlation_plots.Rdata"))
#load(paste0("./Analysis/after_split/Correlations/Selected.Pathways/ICR_score/Selected.pathways_CCSK_Correlation_plots.Rdata"))

TARGET.cancersets = rbind(TARGET.cancersets, c("Pancancer", "Pancancer"))

row.names(TARGET.cancersets) = TARGET.cancersets$cancerType
pancancer_Geneset_cor_table = t(TARGET.cancersets)[-c(1,2),]

pancancer_Geneset_cor_table = rbind(pancancer_Geneset_cor_table,matrix(nrow = nrow(Hallmark_GSEA_cor),ncol=N.sets))

rownames(pancancer_Geneset_cor_table) = rownames(Hallmark_GSEA_cor)

i=1
for (i in 1:N.sets) {
  Cancer = CancerTYPES[i]
  load(paste0("./Analysis/after_split/Correlations/",selected_genes,"/",score,"/",selected_genes, "_",Cancer,"_Correlation_plots.Rdata"))
  if(display_correlations == "only_significant"){
    N.columns = ncol(Hallmark_GSEA_cor)
    N.rows = nrow(Hallmark_GSEA_cor)
    for (i in 1:N.columns){
      for (j in 1:N.rows){
        Hallmark_GSEA_cor[i,j] = ifelse(sele[[1]][i,j] <0.05 | Sel.Path_cor_sign[[1]][i,j] >0.95, Hallmark_GSEA_cor[i,j], 0)
      }
    }
  }
  pancancer_Geneset_cor_table[, Cancer] = as.numeric(Hallmark_GSEA_cor[,score])
}

write.csv(pancancer_Geneset_cor_table,file = "./Analysis/after_split/Correlations/Figure.Supp3.OC_NBL_Pancancer_correlation_with.ICR.csv",quote = FALSE)
# convert to numeric matrix
mode(pancancer_Geneset_cor_table) = "numeric"
#rownames_order = rownames(pancancer_Geneset_cor_table)[order(rowMeans(pancancer_Geneset_cor_table),decreasing = TRUE)]
pancancer_Geneset_cor_table = pancancer_Geneset_cor_table[-which(rownames(pancancer_Geneset_cor_table) == score),]
rownames(pancancer_Geneset_cor_table) = gsub(".*] ", "",rownames(pancancer_Geneset_cor_table))

dev.new()
#Correlation complex heatmap
dir.create(paste0("./Figures/after_split/Correlation_plots/",Gene.Set,"_Correlation_plots") , showWarnings = FALSE)
dir.create(paste0("./Figures/after_split/Correlation_plots/",Gene.Set,"_Correlation_plots/",score),showWarnings = FALSE)
png(paste0("./Figures/after_split/Correlation_plots/",Gene.Set,"_Correlation_plots/",score,"/051.",Gene.Set,"_pearson_",score,"_PanCancer_pediatric.png"),
    res=400,height= 20,width=15,unit="in")

Geneset_cor = pancancer_Geneset_cor_table
cex.before <- par("cex")
par(cex = 0.35)
lims=c(-1,1)
if (length(Geneset_cor[Geneset_cor<0]) == 0) {lims=c(0,1)}
annotation = data.frame (Cancer = colnames(Geneset_cor),color = c("brown" , "#1ABC9C" ,"#4B0082"),stringsAsFactors = FALSE)
#annotation = annotation[corrMatOrder(Geneset_cor,order="FPC"),]
Geneset_cor = as.data.frame(Geneset_cor)
annotation$Cancer = factor(annotation$Cancer, levels = c("Pancancer","OS","NBL_mycn_Namp_High" ))
ha = HeatmapAnnotation(`Cancer` = annotation$Cancer,
                       col = list(`Cancer` = c("NBL_mycn_Namp_High" ="#4B0082" , "OS"= "#1ABC9C" ,"Pancancer"="brown")),annotation_name_gp = gpar(fontsize = 22, fontface = "bold"))

dev.new()
Geneset_cor = as.matrix(Geneset_cor)
HM = Heatmap(Geneset_cor,
             #column_title = paste0("Pearson correlation between ",score," and \n ES of Selected oncogenic pathways signature."),
             column_title_gp = gpar(fontsize = 20, fontface = "bold"), 
             row_title_gp = gpar(fontsize = 10),row_names_gp = gpar(fontsize = 22), # fontface = "bold") ,
             cluster_rows = TRUE ,cluster_columns = FALSE ,
             show_column_names = FALSE,top_annotation = ha,name = "value",col=my.palette,
             # heatmap_legend_param =list(title_gp=gpar(fontsize=10, fontface="bold"),legend_width=unit(8,"cm"),legend_position = "left"),
             #row_names_max_width = unit(5, "cm")
             # theme(legend.position = "none")
             column_order = c("Pancancer","OS","NBL_mycn_Namp_High"),
             row_names_max_width = unit(9, "in")
)
draw(HM, heatmap_legend_side = "left", annotation_legend_side = "left")

par(cex = cex.before)
dev.off()
