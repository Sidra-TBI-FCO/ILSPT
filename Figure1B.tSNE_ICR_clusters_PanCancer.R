#tSNE of ICR CLUSTERS
#Setup environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(master.location)
setwd(paste0(master.location,"/TBI-LAB - Project - Pediatric Pan Cancer TARGET"))
# Install packages and load
source(paste0(toolbox.path,"/R scripts/ipak.function.R"))

required.packages <- c("Rtsne", "ggplot2","dplyr")
ipak(required.packages)  

# Set parameters
Type = "Normalized" # "Raw" or "Normalized"
dataset = "CancerTypes" 
Cancer = "PanCancer"
Gene_set = "HML.classification"

#load data
load(paste0("./Processed_Data/003_PanCancer_pediatric_only_normalized_TCGAbiolinks_TP_filtered.Rdata"))
colnames(filtered.norm.RNAseqData)
filtered.norm.RNAseqData =as.data.frame(filtered.norm.RNAseqData)

#Load The Annotation File

load("./Analysis/after_split/ICR_data_percancer_clustered_408samples/TARGET_percancer_table_cluster.Rdata")


#tSNE
RNASeq.QN.counts = filtered.norm.RNAseqData
RNASeq_QN = t(RNASeq.QN.counts)
RNASeq_QN_LOG2 = log(RNASeq_QN +1, 2)
set.seed(7)

tsne_model_1 = Rtsne(RNASeq_QN_LOG2, check_duplicates=FALSE, pca=TRUE, perplexity=15, theta=0.5, dims=2)
dir.create("./Analysis/after_split/tSNE",showWarnings = FALSE)

#Create the Annotation file ## Add hte percancer clusters file

clustering$HML.ICR.Cluster = factor(clustering$HML.ICR.Cluster, levels = c("ICR High", "ICR Medium", "ICR Low"))

color_table <- tibble(
  ICR_Clusters = c("ICR High", "ICR Medium", "ICR Low"),
  Color = c("red","green","blue"))
d_tsne = as.data.frame(tsne_model_1$Y)

#dev.new()
# Create the tSNE plot 
dir.create("./Figures/after_split/048.tSNE_plots",showWarnings = F, recursive = T)
png("./Figures/after_split/tSNE_plots/048.tSNE_ICR_percancer_clustering.png", height = 4, width = 6,
    units = "in", res = 600)
plot= ggplot(d_tsne, aes(x=V1, y=V2,color=clustering$HML.ICR.Cluster)) +
  geom_point(aes(x=V1, y=V2),size=0.7) +
  guides(colour=guide_legend(override.aes=list(size=6),title = "ICR Clusters")) +
  #scale_colour_brewer(palette = "Set2")+
  scale_color_manual(values=c(color_table$Color))+
  xlab("tSNE dimension 1") + ylab("tSNE dimension 2") +
  #ggtitle("Perplexity=15") +
  theme_light(base_size=20) +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        #legend.position 
  )
plot(plot)
dev.off()
