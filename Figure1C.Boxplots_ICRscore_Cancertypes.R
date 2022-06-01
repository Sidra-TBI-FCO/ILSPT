#Calculation of ICR score
# Install packages and load
rm(list=ls())
load("~/R.Config.Rdata")
setwd(master.location)
setwd(paste0(master.location,"/TBI-LAB - Project - Pediatric Pan Cancer TARGET"))
source(paste0(toolbox.path,"/R scripts/ipak.function.R"))
required.packages = c("stringr", "ggplot2", "ggpubr")
ipak(required.packages)

## Set Parameters
Cancer = "PanCancer"
signature = c("IFNG","IRF1", "STAT1","IL12B","TBX21","GZMB","GNLY","PRF1","GZMH", "GZMA","CXCL10","CXCL9","CCL5","CD8A","CD8B","PDCD1","CTLA4","FOXP3","IDO1","CD274") #ICR 

classification_k = "Clusters"
score = "ICRscore"
basis_ordering = "Mean ICR"
# Load data
load("./Analysis/after_split/annotation_with_NBL_subgroups.Rdata")

load(paste0("./Processed_Data/003_PanCancer_pediatric_only_normalized_TCGAbiolinks_TP_filtered.Rdata"))
filtered.norm.RNAseqData = filtered.norm.RNAseqData[which(rownames(filtered.norm.RNAseqData) %in% signature),]
filtered.norm.RNAseqData = filtered.norm.RNAseqData[,which(colnames(filtered.norm.RNAseqData) %in% annotation_all$Sample)]
subset_RNASeq = t(filtered.norm.RNAseqData)
selected_subset_RNAseq_log2 = log(subset_RNASeq +1, 2)
annotation_all[,score] = rowMeans(selected_subset_RNAseq_log2) 

dir.create("./Figures/after_split", showWarnings = FALSE)
#dir.create(paste0("./Figures/after_split/Boxplots/"), showWarnings = FALSE)
annotation = annotation_all

annotation$Cancer = factor(annotation$Cancer,levels = c("CCSK","OS","RT","WT", "NBL_mycn_Amp","NBL_mycn_Namp_Intermed.Low","NBL_mycn_Namp_High"))

annotation = annotation[!is.na(annotation[,score]),]

if(basis_ordering == "Mean ICR"){
  mean_ICR = aggregate(annotation[,score], list(annotation$Cancer), mean)
  mean_ICR = mean_ICR[order(mean_ICR$x, decreasing = FALSE),]
  ICR_order = mean_ICR$Group.1
}

annotation$Cancer = factor(annotation$Cancer,levels = ICR_order)


dev.new()
#PLOTTING
png(paste0("./Figures/after_split/Boxplots/037.CancerTypes.plus.NBL.subgroups_",score,"_boxplot.png"), res=400,height=6,width=7,unit="in")
plot = ggplot(annotation, aes(x= Cancer, y= ICRscore, fill = Cancer)) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(size = 0.1, width = 0.15) +
  theme_bw() +
  scale_fill_manual(values = c("CCSK" = "#BB8FCE","OS"="#1ABC9C","RT"="#FA8072","WT"="#2980B9","NBL_mycn_Amp"="#FF4500","NBL_mycn_Namp_Intermed.Low"="#D78414","NBL_mycn_Namp_High"="#4b0082")) +
  scale_y_continuous("ICR score") +
  xlab("Cancer") +
  theme(
    #axis.text.x = element_text(color = "black", size = 10),
        axis.text.x = element_text(angle=50, hjust=1),
        axis.text.y = element_text(color = "black", size = 10),
        axis.title.x = element_text(color = "black", size = 15),
        axis.title.y = element_text(color = "black", size = 15))
plot(plot)
dev.off()
