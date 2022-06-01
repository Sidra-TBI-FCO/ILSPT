
# Script for generating the Boxplot of ICRscore across MYCN groups of NBL

#Setup microenvironment
rm(list = ls())
setwd("~/Sidra Medicine - Research Division/TBI-LAB - Project - Pediatric Pan Cancer TARGET/")

# Install packages and load
source("~/Sidra Medicine - Research Division/TBI-LAB - Project - Pediatric Pan Cancer TARGET/tools/ipak.function.R")
required.packages = c("stringr", "ggplot2", "ggpubr")
ipak(required.packages)

## Set Parameters
ICR_classification_k = "MYCN.status"
basis_ordering = "Mean ICR"
subset = "all"                                                                                                         
Cancer  = c("NBL_mycn_Namp_High", "NBL_mycn_Namp_Intermed.Low" , "NBL_mycn_Amp")

# Load data
load("./Analysis/after_split/ICR_data_Pancancer/TARGET_Pancancer_table_cluster.RData")
clinical_NBL= read.csv("./Data/Clinical_Data/046.TARGET-AllCancers_clinical.plus.NBL.csv",stringsAsFactors = FALSE) 

clinical_NBL = clinical_NBL[which(clinical_NBL$Cancer %in% Cancer),]
clustering = clustering[which(rownames(clustering) %in% clinical_NBL$submitter_id),]

dir.create(paste0("./Figures/after_split/Boxplots/MYCN.status_NBL151_boxplot_across_cancers/"), showWarnings = FALSE)
clinical_NBL$ICRscore = clustering$ICRscore[match(clinical_NBL$submitter_id,rownames(clustering))]
clinical_NBL = clinical_NBL[!is.na(clinical_NBL$ICRscore),]

known_levels <- c("NBL_mycn_Namp_Intermed.Low","NBL_mycn_Amp","NBL_mycn_Namp_High")
my_order <- order(factor(clinical_NBL$Cancer, levels = known_levels, ordered=TRUE))
clinical_NBL = clinical_NBL[my_order, ]


clinical_NBL$Cancer <- factor(clinical_NBL$Cancer,levels = c("NBL_mycn_Namp_Intermed.Low","NBL_mycn_Amp","NBL_mycn_Namp_High"))


dev.new()
png(paste0("./Figures/after_split/Boxplots/MYCN.status_NBL151_boxplot_across_cancers/018.Violin_plots_colour_by_ICR_cluster_across_NBL_mycn_groups.png"), res=400,height=6,width=9,unit="in")
plot = ggplot(clinical_NBL, aes(x = Cancer, y = ICRscore, fill = Cancer)) +
  geom_violin(outlier.shape=NA) +
  geom_boxplot(width = 0.3) +
  geom_jitter(size = 0.1, width = 0.15) +
  theme_bw() +
  scale_fill_manual(values = c("#D78414","#FF4500","#4b0082")) +
  scale_y_continuous("ICRscore") +
  xlab("Cancer") +
  stat_compare_means(method = "t.test", comparisons = list(c("NBL_mycn_Namp_High", "NBL_mycn_Namp_Intermed.Low"),
                                                           c("NBL_mycn_Namp_High","NBL_mycn_Amp"),
                                                           c("NBL_mycn_Namp_Intermed.Low","NBL_mycn_Amp"))) +
  theme(axis.text.x = element_text(color = "black", size = 10),
        axis.text.y = element_text(color = "black", size = 10),
        axis.title.x = element_text(color = "black", size = 10),
        axis.title.y = element_text(color = "black", size = 16))

print(plot)
dev.off()
