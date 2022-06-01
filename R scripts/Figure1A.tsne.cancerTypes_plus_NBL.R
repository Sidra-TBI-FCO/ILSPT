# Script for generating tSNE plot annotated by cancer types

#Setup environment
rm(list = ls())

# Install packages and load
required.packages <- c("Rtsne", "ggplot2","dplyr")
library(required.packages)  

# Set parameters
Type = "Normalized" # "Raw" or "Normalized"
dataset = "Cancertypes" 
Cancer = "PanCancer"
Gene_set = "filtered.normalized"

#load data
load(paste0("./Processed_Data/003_PanCancer_pediatric_only_normalized_TCGAbiolinks_TP_filtered.Rdata"))
filtered.norm.RNAseqData = as.data.frame(filtered.norm.RNAseqData)

# Load the annotation file
load("./Analysis/after_split/annotation_with_NBL_subgroups.Rdata")
annotation = annotation_all

#tSNE
RNASeq.QN.counts = filtered.norm.RNAseqData
RNASeq_QN = t(RNASeq.QN.counts)
RNASeq_QN_LOG2 = log(RNASeq_QN +1, 2)
set.seed(7)

tsne_model_1 = Rtsne(RNASeq_QN_LOG2, check_duplicates=FALSE, pca=TRUE, perplexity=15, theta=0.5, dims=2)

# Create the output directory
dir.create("./Analysis/after_split/tSNE",showWarnings = FALSE)

annotation$Cancer = factor(annotation$Cancer, levels = c("CCSK","OS","RT","WT", "NBL_mycn_Amp","NBL_mycn_Namp_Intermed.Low","NBL_mycn_Namp_High"))
color_table <- tibble(
  Cancer = c("CCSK","OS","RT","WT", "NBL_mycn_Amp","NBL_mycn_Namp_Intermed.Low","NBL_mycn_Namp_High"),
  Color = c("#BB8FCE","#1ABC9C","#FA8072","#2980B9","#FF4500","#D78414","#4b0082"))
d_tsne = as.data.frame(tsne_model_1$Y)

# Create the tSNE plot 
dir.create("./Figures/after_split/tSNE_plots",showWarnings = F, recursive = T)
png("./Figures/after_split/tSNE_plots/049.tSNE_CancerTypes.PLUS.nbl.subtypes.png", height = 4, width = 8,
    units = "in", res = 600)

plot= ggplot(d_tsne, aes(x=V1, y=V2,color=annotation$Cancer)) +
  geom_point(aes(x=V1, y=V2),size=.7) +
  guides(colour=guide_legend(override.aes=list(size=6),title = "Cancer types")) +
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
