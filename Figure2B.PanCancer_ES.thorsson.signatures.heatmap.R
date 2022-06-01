#Enrichment scores HEATMAP (5 thorsson clusters)
#Setup environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(master.location)
setwd(paste0(master.location,"/TBI-LAB - Project - Pediatric Pan Cancer TARGET"))
# Install packages and load
source(paste0(toolbox.path,"/R scripts/ipak.function.R"))
required.packages = c("ComplexHeatmap","factoextra")
ipak(required.packages)

# Set parameters
Cancer = "PanCancer"
Gene.set = "thorss_list_5_clusters"

# load data
load(paste0("./Analysis/after_split/Signature_Enrichment/GSEA_PanCancer_pediatric_onlythorss_list_5_clusters.Rdata"))

load(paste0("./Analysis/after_split/annotation_with_NBL_subgroups.Rdata"))
ES = ES[,which(colnames(ES) %in% annotation_all$Sample)]

# Create directories
dir.create("./Figures", showWarnings = FALSE)
dir.create(paste0("./Figures/after_split/Heatmaps/", Gene.set), showWarnings = FALSE)

# z-score Expression.matrix
ESz = ES

for(j in 1: nrow(ESz))  {
  ESz[j,] = (ES[j,]-mean(ES[j,]))/sd(ES[j,])
}

rownames(ESz)[rownames(ESz) == "CHANG_CORE_SERUM_RESPONSE_UP"] <- "Wound healing"
rownames(ESz)[rownames(ESz) == "Module3_IFN_score"] <- "IFN_G"
rownames(ESz)[rownames(ESz) == "CSF1_response"] <- "Macrophage"
rownames(ESz)[rownames(ESz) == "TGFB_score_21050467"] <- "TGF_B"
rownames(ESz)[rownames(ESz) == "LIexpression_score"] <- "Lymphocyte"

annotation_all = data.frame(Sample = annotation_all$Sample,
                        Cancer = annotation_all$Cancer)
rownames(annotation_all) = annotation_all$Sample
ESz = ESz[,which(colnames(ESz) %in% annotation_all$Sample)]

annotation_all = annotation_all[colnames(ESz),]

row_ha = rowAnnotation(
  Cancer = annotation_all$Cancer,
  #KM_5 = annotation2$cluster,
  col = list(`Cancer` = c("CCSK" = "#BB8FCE","OS"="#1ABC9C","RT"="#FA8072","WT"="#2980B9","NBL_mycn_Amp"="#FF4500","NBL_mycn_Namp_Intermed.Low"="#D78414","NBL_mycn_Namp_High"="#4b0082")),
  block = anno_block(gp = gpar(fill = c("#DF536B","#7ECD60","#4A93DF","#6DDFE3","#BC3DB5","#FFA500"),labels = c("S1", "S2", "S3","S4","S5","S6"),
                               labels_gp = gpar(col = "white", fontsize = 10))), 
  width = unit(12, "mm"))


ESz = t(ESz)
colnames(ESz)

col_order <- c("IFN_G", "TGF_B", "Macrophage",   #reorder the column names as in thorsson paper 
               "Lymphocyte", "Wound healing")
ESz <- ESz[, col_order]

ESz1 = scale(ESz)
######
#To select the optima number of clusters 
#set.seed(1111)
#fviz_nbclust(ESz, kmeans, method = "wss")

#fviz_nbclust(ESz, kmeans, method = "wss") +
 # geom_vline(xintercept = 4, linetype = 2)+
  #labs(subtitle = "Elbow method")
#############################################

png(filename = paste0("./Figures/after_split/Heatmaps/",Gene.set,"/030.Complexheatmap_", Cancer ,"_", Gene.set, "_km6_10000repeats_heatmap_Euclidean_pediatric_only_plus_NBL_groups.V3.png"), res = 1500,
    width = 8, height = 6, units = "in")

HM = Heatmap(ESz, cluster_rows = TRUE ,cluster_columns = FALSE , row_names_max_width = unit(5, "in"),
             row_km=6,row_km_repeats = 10000 ,cluster_row_slices = FALSE,clustering_distance_rows = "euclidean",
             show_column_names = TRUE,row_title = "S%s", show_row_names = FALSE, left_annotation = row_ha,name = "Enrichment\n z score",
             #row_names_gp = gpar(fontsize = 0.4)
             )
print(HM)
dev.off()

# to rename the clusters 
library(circlize)
 
r.dend <- row_dend(HM)  #Extract row dendogram
#set.seed(104)
rcl.list <- row_order(HM)  #Extract clusters (output is a list)

cluster_1_samples = rownames(ESz)[rcl.list[["1"]]]
cluster_2_samples = rownames(ESz)[rcl.list[["2"]]]
cluster_3_samples = rownames(ESz)[rcl.list[["3"]]]
cluster_4_samples = rownames(ESz)[rcl.list[["4"]]]
cluster_5_samples = rownames(ESz)[rcl.list[["5"]]]
cluster_6_samples = rownames(ESz)[rcl.list[["6"]]]

annotation_all$cluster = NA
annotation_all$cluster[which(annotation_all$Sample %in% cluster_1_samples)] = "1"
annotation_all$cluster[which(annotation_all$Sample %in% cluster_2_samples)] = "2"
annotation_all$cluster[which(annotation_all$Sample %in% cluster_3_samples)] = "3"
annotation_all$cluster[which(annotation_all$Sample %in% cluster_4_samples)] = "4"
annotation_all$cluster[which(annotation_all$Sample %in% cluster_5_samples)] = "5"
annotation_all$cluster[which(annotation_all$Sample %in% cluster_6_samples)] = "6"

table(annotation_all$Cancer, annotation_all$cluster)

annotation_all$cluster[annotation_all$cluster == "1"] <-"S1"
annotation_all$cluster[annotation_all$cluster == "2"] <-"S2"
annotation_all$cluster[annotation_all$cluster == "3"] <-"S3"
annotation_all$cluster[annotation_all$cluster == "4"] <-"S4"
annotation_all$cluster[annotation_all$cluster == "5"] <-"S5"
annotation_all$cluster[annotation_all$cluster == "6"] <-"S6"

save(annotation_all,file ="./Analysis/after_split/Signature_Enrichment/Thorsson_clustering/030.Thorsson_samples_cancers_10000repeats_HM_clusters_km6_pediatric_only_V3.Rdata")
