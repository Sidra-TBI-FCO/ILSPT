#Script for generating the RNASeq complex heatmap 

# Setup environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(master.location)
setwd(paste0(master.location,"/TBI-LAB - Project - Pediatric Pan Cancer TARGET"))
# Install packages and load
source(paste0(toolbox.path,"/R scripts/ipak.function.R"))
required.packages = c("ComplexHeatmap")
ipak(required.packages)

# Set parameters
groups = "HML.ICR.Cluster" # "k4" or "HL_k4" or "k2" or "HL_k2"
Cancer = "PanCancer"

# load data
load("./Analysis/after_split/Signature_Enrichment/Thorsson_clustering/030.Thorsson_samples_cancers_10000repeats_HM_clusters_km6_pediatric_only_V3.Rdata")
load(paste0("./Processed_Data/003_PanCancer_pediatric_only_normalized_TCGAbiolinks_TP_filtered.Rdata"))
load("~/Sidra Medicine - Research Division/TBI-LAB - General/bioinformatics tools/ICR genes/ICR_genes.RData")
load(paste0("./Analysis/after_split/ICR_data_percancer_clustered_408samples/TARGET_percancer_table_cluster.Rdata"))

clustering$Cancer = annotation_all$Cancer[match(rownames(clustering),annotation_all$Sample)]

# Create directories
#dir.create("./Analysis/after_split/", showWarnings = FALSE)
#dir.create(paste0("./Analysis/after_split/ICR_data",Cancer), showWarnings = FALSE)

# subset expression matix for ICR
ICR_subset_RNAseq = filtered.norm.RNAseqData[which(row.names(filtered.norm.RNAseqData) %in% ICR_genes), ] 
ICR_subset_RNAseq_log2 = log(ICR_subset_RNAseq +1, 2)

clustering = clustering[order(clustering$ICRscore),]
sample_order = rownames(clustering)

Expression.matrix = ICR_subset_RNAseq_log2[,sample_order]

# z-score Expression.matrix
Expression.matrix.z = Expression.matrix
for(j in 1: nrow(Expression.matrix.z))  {
  Expression.matrix.z[j,] = (Expression.matrix[j,]-mean(Expression.matrix[j,]))/sd(Expression.matrix[j,]) # z-score the enrichment matrix
}

clustering$HML.ICR.Cluster = factor(clustering$HML.ICR.Cluster, levels = c("ICR High", "ICR Medium", "ICR Low"))


annotation = data.frame(Sample = rownames(clustering),
                        ICR_Cluster = clustering$HML.ICR.Cluster,
                        Cancer = clustering$Cancer)

rownames(annotation) = annotation$Sample

annotation = annotation[colnames(Expression.matrix.z),]

ha = HeatmapAnnotation(ICR_cluster = annotation$ICR_Cluster,
                       Cancer = annotation$Cancer,
                       col =  list(`ICR_cluster`= c("ICR High" = "red", "ICR Medium" = "green", "ICR Low" = "blue"),
                                   `Cancer`      = c("CCSK" = "#BB8FCE", "NBL_mycn_Namp_High" = "#4b0082","NBL_mycn_Namp_Intermed.Low" = "#D78414", "NBL_mycn_Amp"="#FF4500", "WT" = "#2980B9" , "OS"= "#1ABC9C" , "RT"="#FA8072" )))


png(filename = paste0("./Figures/after_split/Heatmaps/ICR_heatmaps/007.",Cancer,"_ICR_ComplexHeatmap_", groups, "_annotated_by_CancerTypes.png"), res = 600,
    width = 9, height = 6, units = "in")

dev.new()
Heatmap(Expression.matrix.z, cluster_rows = TRUE ,cluster_columns = TRUE , row_names_max_width = unit(5, "in"),
        show_column_names = FALSE, top_annotation =ha, name = "Expression\n z score"
)

dev.off()
