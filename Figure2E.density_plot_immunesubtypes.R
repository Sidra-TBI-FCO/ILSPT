# Script for generating the density plot of immune signatures

# Setup environemnt
rm(list = ls())
load("~/R.Config.Rdata")

setwd(paste0("~/Sidra Medicine - Research Division/TBI-LAB - Project - Pediatric Pan Cancer TARGET/"))

source(paste0("../TBI-LAB - General/Bioinformatics tools/R scripts/ipak.function.R"))

required.packages = c("colorspace", "ggplot2","dplyr","ComplexHeatmap","factoextra","reshape2")
ipak(required.packages)

#Parameters 
Cancer = "PanCancer"
Gene.set = "Thorsson"

#Loading the required data
load("./Analysis/after_split/Signature_Enrichment/GSEA_PanCancer_pediatric_onlyThorsson.Rdata")
load("./Analysis/after_split/Signature_Enrichment/Thorsson_clustering/030.Thorsson_samples_cancers_10000repeats_HM_clusters_km6_pediatric_only_V3.Rdata")
ESz = ESz[,which(colnames(ESz) %in% annotation_all$Sample)]
ES = ES[,which(colnames(ES) %in% annotation_all$Sample)]

rownames(ESz)[rownames(ESz) == "Wolf.CHANG_CORE_SERUM_RESPONSE_UP"] <- "Wound_healing"
rownames(ESz)[rownames(ESz) == "Wolf.Module3_IFN_score"] <- "IFN_G"

rownames(ESz)[rownames(ESz) == "Wolf.CSF1_response"] <- "Macrophage"
rownames(ESz)[rownames(ESz) == "Wolf.TGFB_score_21050467"] <- "TGF_B"
rownames(ESz)[rownames(ESz) == "Wolf.LIexpression_score"] <- "Lymphocyte"
rownames(ESz)[rownames(ESz) == "Bindea.Th1 cells"] <- "Th1_cells"
rownames(ESz)[rownames(ESz) == "Bindea.Th2 cells"] <- "Th2_cells"
rownames(ESz)[rownames(ESz) == "Bindea.Cytotoxic cells"] <- "Cytotoxic_cells"
rownames(ESz)[rownames(ESz) == "Wolf.Module11_Prolif_score"] <- "Proliferation"
rownames(ESz)[rownames(ESz) == "Bindea.Th17 cells"] <- "Th17"
rownames(ESz)[rownames(ESz) == "Bindea.Treg cells"] <- "Treg"
rownames(ESz)[rownames(ESz) == "Bindea.NK cells"] <- "NK_cells"
rownames(ESz)[rownames(ESz) == "Bindea.NK CD56bright cells"] <- "NK_cells_CD56bright"
rownames(ESz)[rownames(ESz) == "Bindea.NK CD56dim cells"] <- "NK_cells_CD56dim"
rownames(ESz)[rownames(ESz) == "ICR.ICR_SCORE"] <- "ICR_SCORE"

# Adding the HLA 1,2,3 data 
Gene.set = c("HHLA1","HHLA2","HHLA3","HLA-DPB2","HLA-DRB5","HLA-DRB6","HLA-F","HLA-F-AS1","HLA-H","ICR_SCORE")
load(paste0("./Processed_Data/003_PanCancer_pediatric_only_normalized_TCGAbiolinks_TP_filtered.Rdata"))
filtered.norm.RNAseqData = filtered.norm.RNAseqData[,which(colnames(filtered.norm.RNAseqData) %in% annotation_all$Sample)]
load( "./Analysis/after_split/Signature_Enrichment/Thorsson_clustering/030.Thorsson_samples_cancers_10000repeats_HM_clusters_km6_pediatric_only_V3.Rdata")
filtered.norm.RNAseqData = filtered.norm.RNAseqData[which(rownames(filtered.norm.RNAseqData) %in% Gene.set),]
filtered.norm.RNAseqData = log(filtered.norm.RNAseqData +1, 2)
filtered.norm.RNAseqData = t(filtered.norm.RNAseqData)
filtered.norm.RNAseqData = as.data.frame(filtered.norm.RNAseqData)
ESz = as.data.frame(t(ESz))

ESz$HLA1 = filtered.norm.RNAseqData$HHLA1[match(rownames(ESz),rownames(filtered.norm.RNAseqData))]
ESz$HLA2 = filtered.norm.RNAseqData$HHLA2[match(rownames(ESz),rownames(filtered.norm.RNAseqData))]
ESz$HLA3 = filtered.norm.RNAseqData$HHLA3[match(rownames(ESz),rownames(filtered.norm.RNAseqData))]

ESz = t(ESz)

signatures = c("Wound_healing","IFN_G","Macrophage","TGF_B","Lymphocyte","Th1_cells","Th2_cells","Th17","Cytotoxic_cells","Treg","HLA1","HLA2","ICR_SCORE")


ESz = ESz[which(rownames(ESz) %in% signatures),]

ESz = t(ESz)

ESz = as.data.frame(ESz)
colMeans(ESz)

#ESz = scale(ESz,scale = TRUE,center = TRUE )

ESz = ESz[which(rownames(ESz) %in% annotation_all$Sample),]
ESz$cluster = annotation_all$cluster[match(annotation_all$Sample,rownames(ESz))]

#########
dir.create("./Figures/after_split/Density_plots",showWarnings = FALSE)
#png(filename = paste0("./Figures/after_split/Density_plots/Density_plot_", Cancer ,"_", Gene.set, "_Thorsson_5_immun_subtypes.TEST.png"), res = 1500,
 #   width = 15, height = 12, units = "in")

ESz$sampleID =rownames(ESz)
ESz.melted = melt(ESz,id.vars = c("sampleID","cluster"))
ESz.melted = ESz.melted %>% group_by(cluster,variable) %>%
  mutate(median = median(value))
unique(ESz.melted$variable)

ESz.melted$variable = factor(ESz.melted$variable,levels = c("IFN_G","TGF_B","Macrophage","Lymphocyte","Wound_healing","Th1_cells","Th2_cells","Th17","Treg","Cytotoxic_cells","HLA1","HLA2","ICR_SCORE"))
#ESz.melted$variable = factor(ESz.melted$variable,levels = c("IFN_G","TGF_B","Macrophage","Lymphocyte","Wound_healing","Th1_cells","Th2_cells","Th17","Treg","Cytotoxic_cells","HLA1","HLA2","Proliferation","NK_cells","NK_cells_CD56bright","NK_cells_CD56dim"))

dev.new()
median(ESz.melted$value[which(ESz.melted$variable == "Lymphocyte")])

png(filename = paste0("./Figures/after_split/Density_plots/026.Density_plot_", Cancer ,"_", Gene.set, "_km6_12.other_signatures.with.ICR.png"), res = 600,
    width = 10, height = 4.5, units = "in")

ggplot(data = ESz.melted,aes(x=value)) +
  geom_density(aes(fill=median), stat = "density",) +
  scale_fill_gradient2(low ="blue", mid ="white",high ="red") +
  geom_vline(aes(xintercept= median(value)),
             color="black", linetype="dashed", size=0.5) +
  facet_grid(cluster ~ variable,scales = "free")+
  ylim(c(0,1.4)) + xlim(c(-3,6)) +
theme(axis.text.x = element_text(size = 8,color = "black"))+
theme(strip.background = element_blank(), strip.text = element_blank(),plot.background =element_blank(),panel.border = element_blank(),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),
    axis.text.y = element_blank(),axis.ticks.y = element_blank(),
  # axis.text.x = element_text(),axis.ticks.x = element_blank() # Without titles at all 
)
  dev.off()
  
