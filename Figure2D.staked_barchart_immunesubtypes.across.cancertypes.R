#Script for generating Stacked barchart for the clusters across cancer types 
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
load( "./Analysis/after_split/Signature_Enrichment/Thorsson_clustering/030.Thorsson_samples_cancers_10000repeats_HM_clusters_km6_pediatric_only_V3.Rdata")

annotation_all$Cancer = as.character(annotation_all$Cancer)

annotation_all$Cancer[which(annotation_all$Cancer == "NBL_mycn_NA_High")] = "NBL_mycn_Namp_High"
annotation_all$Cancer[which(annotation_all$Cancer == "NBL_mycn_NA_Intermed.low")] = "NBL_mycn_Namp_Intermed.low"

DF <- annotation_all %>%
  group_by(Cancer, cluster) %>%
  summarise(count=n()) %>%
  mutate(perc=count/sum(count))


dir.create(paste0("./Figures/after_split/stacked_Bar_chart"),showWarnings = FALSE)
dir.create(paste0("./Figures/after_split/stacked_Bar_chart/",Gene.set),showWarnings = FALSE)
png(filename = paste0("./Figures/after_split/stacked_Bar_chart/",Gene.set,"/027.", Cancer ,"_", Gene.set, "_HM_clustering_KM6_10000repeats_V3.png"), res = 400,
    width = 5, height = 5, units = "in")
plot = ggplot(DF, aes(x = Cancer, y =perc*100, fill = cluster)) + geom_bar(stat="identity") +
  labs(x = "Cancer", y = "Percentage", fill = "cluster", face = "bold") +
  theme(panel.grid = element_line(linetype = "solid", colour = "white"), 
        panel.background = element_rect(fill = "white", colour = "grey"),
        axis.text.x = element_text(size = 10, colour = "black", angle = 35,
                                   vjust = 1, hjust = 1),
        axis.text.y = element_text(size = 19, colour = "black"),
        axis.title = element_text(size = 19, colour = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 19, colour = "black")) +
  scale_fill_manual(values= c("S1"="#DF536B","S2"="#7ECD60","S3"="#4A93DF","S4"="#6DDFE3","S5"="#BC3DB5","S6"="#FFA500"))
plot(plot)
dev.off()

#c("CCSK" = "#BB8FCE", "NBL" = "#FFC300","OS"= "#1ABC9C","RT"="#FA8072","WT"= "#2980B9")