# Script for generating the stacked barchart for the cancer types  within the 6 immune subtypes 

#Setup environment
rm(list = ls())
load("~/R.Config.Rdata")

setwd(master.location)
setwd(paste0(master.location,"/TBI-LAB - Project - Pediatric Pan Cancer TARGET"))
source(paste0(toolbox.path,"/R scripts/ipak.function.R"))

required.packages = c("stringr", "ggplot2", "ggpubr", "dplyr")
ipak(required.packages)

# Set parameters
Gene.set = "thorss_list_5_clusters"
Cancer = "PanCancer"

# Load data
load("./Analysis/after_split/Signature_Enrichment/Thorsson_clustering/030.Thorsson_samples_cancers_10000repeats_HM_clusters_km6_pediatric_only_V3.Rdata")
#annotation_all$Cancer = as.character(annotation_all$Cancer)

#annotation_all$Cancer[which(annotation_all$Cancer == "NBL_mycn_NA_High")] = "NBL_mycn_Namp_High"
#annotation_all$Cancer[which(annotation_all$Cancer == "NBL_mycn_Namp_Intermed.low")] = "NBL_mycn_Namp_Intermed.Low"

#save(annotation_all,file = "./Analysis/after_split/Signature_Enrichment/Thorsson_clustering/030.Thorsson_samples_cancers_10000repeats_HM_clusters_km6_pediatric_only_V3.Rdata")

DF <- annotation_all %>%
  group_by(cluster, Cancer) %>%
  summarise(count=n()) %>%
  mutate(perc=count/sum(count))


#dev.new()
dir.create(paste0("./Figures/after_split/stacked_Bar_chart"),showWarnings = FALSE)
dir.create(paste0("./Figures/after_split/stacked_Bar_chart/",Gene.set),showWarnings = FALSE)
png(filename = paste0("./Figures/after_split/stacked_Bar_chart/",Gene.set,"/028.", Cancer ,"_", Gene.set, "_cluster_HM_clustering_KM6_10000repeats_V3.png"), res = 400,
    width = 9, height = 5, units = "in")
plot = ggplot(DF, aes(x = cluster, y =perc*100, fill = Cancer)) + geom_bar(stat="identity") +
  labs(x = "cluster", y = "Percentage", fill = "Cancer", face = "bold") +
  theme(panel.grid = element_line(linetype = "solid", colour = "white"), 
        panel.background = element_rect(fill = "white", colour = "grey"),
        axis.text.x = element_text(size = 10, colour = "black", angle = 35,
                                   vjust = 1, hjust = 1),
        axis.text.y = element_text(size = 19, colour = "black"),
        axis.title = element_text(size = 19, colour = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 19, colour = "black")) +
  scale_fill_manual(values= c("CCSK" = "#BB8FCE","OS"="#1ABC9C","RT"="#FA8072","WT"="#2980B9","NBL_mycn_Amp"="#FF4500","NBL_mycn_Namp_Intermed.Low"="#D78414","NBL_mycn_Namp_High"="#4b0082"))
plot(plot)
dev.off()

#  "CCSK" = "#BB8FCE","OS"="#1ABC9C","RT"="#FA8072","WT"="#2980B9","NBL_mycn_Amp"="#FF4500","NBL_mycn_NA_Intermed.low"="#D78414","NBL_mycn_NA_High"="#4b0082"
#c("CCSK" = "#BB8FCE", "NBL" = "#FFC300","OS"= "#1ABC9C","RT"="#FA8072","WT"= "#2980B9")
