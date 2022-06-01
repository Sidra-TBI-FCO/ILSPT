# Script for generating boxplots for immune checkpoints expression across immune subtypes

# Setup environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(master.location)
setwd(paste0(master.location,"/TBI-LAB - Project - Pediatric Pan Cancer TARGET"))
# Install packages and load
source(paste0(toolbox.path,"/R scripts/ipak.function.R"))
required.packages = c("stringr", "ggplot2", "ggpubr")
ipak(required.packages)

# Set parameters
Cancer = "Pancancer"
Gene.set = "ImmuneCheckPoints"
test = "t.test"
cluster = c("S1","S2","S3","S4","S5","S6")
# Load data
Checkpoints = read.csv("~/Sidra Medicine - Research Division/TBI-LAB - General/bioinformatics tools/Checkpoint genes/checkpoint.genes.csv",stringsAsFactors = FALSE)
load(paste0("./Processed_Data/003_PanCancer_pediatric_only_normalized_TCGAbiolinks_TP_filtered.Rdata"))
load( "./Analysis/after_split/Signature_Enrichment/Thorsson_clustering/Thorsson_samples_cancers_10000repeats_HM_clusters_km6_pediatric_only_V2.Rdata")
#rownames(ES) = gsub("\\/", "_", rownames(ES))

filtered.norm.RNAseqData = filtered.norm.RNAseqData[which(rownames(filtered.norm.RNAseqData) %in% Checkpoints$HUGO.GENE.SYM),]

filtered.norm.RNAseqData = log(filtered.norm.RNAseqData +1, 2)
#annotation = annotation[which(annotation$Cancer %in% Cancer),]
annotation = annotation[which(annotation$cluster %in% cluster),]
filtered.norm.RNAseqData = filtered.norm.RNAseqData[,which(colnames(filtered.norm.RNAseqData) %in% annotation$Sample)]

# Analysis
data = data.frame(Sample_ID = annotation$Sample,
                  cluster = annotation$cluster,
                  Pathway_ES = NA)


data$cluster = factor(data$cluster, levels = c("S1","S2","S3","S4","S5","S6"))


dir.create("./Figures/after_split/Boxplots/thorss_clusters_10000_km6", showWarnings = FALSE)

dir.create(paste0("./Figures/after_split/Boxplots/thorss_clusters_10000_km6/",Gene.set), showWarnings = FALSE)
dir.create(paste0("./Figures/after_split/Boxplots/thorss_clusters_10000_km6/",Gene.set,"/",test,"_BoxPlot/"), showWarnings = FALSE)
dev.new()
#plotting
i = 11 
for (i in 1:nrow(filtered.norm.RNAseqData)){
  pathway = rownames(filtered.norm.RNAseqData)[i]
  data$Pathway_ES = filtered.norm.RNAseqData[pathway,][match(data$Sample_ID, colnames(filtered.norm.RNAseqData))]
  
  plot = ggplot(data, aes(x = cluster, y = Pathway_ES, fill = cluster)) +
    geom_violin(outlier.shape=NA) +
    geom_boxplot(width = 0.3) +
    geom_jitter(size = 0.1, width = 0.15) +
    theme_bw() +
    scale_fill_manual(values = c("#DF536B","#7ECD60","#4A93DF","#6DDFE3","#BC3DB5","#FFA500")) +
    ylab(paste0(pathway, " log2 expression")) +
    xlab("") +
    stat_compare_means(method = "t.test", comparisons = list(c("S2","S4"),
                                                            c("S2","S6"))) +
    theme(axis.text.x = element_text(color = "black", size = 10),
          axis.text.y = element_text(color = "black", size = 10),
          axis.title.x = element_text(color = "black", size = 10),
          axis.title.y = element_text(color = "black", size = 16))
  
  png(file = paste0("./Figures/after_split/Boxplots/thorss_clusters_10000_km6/",Gene.set,"/",test,"_BoxPlot/", Cancer, "_", pathway, "_by_clusters.png"), 
      res = 600, units = "in", width = 5, height = 4)
  plot(plot)
  dev.off()
}
