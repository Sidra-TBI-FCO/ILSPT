# Script for performing Spearman correlation heatmap of ICR Per Cancer

# Setup environment
rm(list = ls())
load("~/R.Config.Rdata")

setwd(master.location)
setwd(paste0(master.location,"/TBI-LAB - Project - Pediatric Pan Cancer TARGET"))

source(paste0(toolbox.path,"/R scripts/ipak.function.R"))
required.packages <- c("corrplot", "stringr")
ipak(required.packages)                                                                                                  # Install and load required packages

# Set Parameters
Cancer = "WT"       ##  	             NBL_mycn_Amp              NBL_mycn_Namp_Intermed.Low                            # Specify the cancertypes that you want to download or process, c("...","...") or "ALL"
colpattern = colorRampPalette(c("#152B7E", "white", "#1B7E09"))(n = 444)
selected_genes = "ICR"                                                                                                   # Specify which genes will be correlated
test = "pearson"

# Load data

#load(paste0("./Processed_Data/TARGET_Pan_cancer/",Cancer,"_PanNormalized_matrix.Rdata"))
#for PanCancer 
load(paste0("./Processed_Data/003_PanCancer_pediatric_only_normalized_TCGAbiolinks_TP_filtered.Rdata"))
load("./Analysis/after_split/Signature_Enrichment/Thorsson_clustering/030.Thorsson_samples_cancers_10000repeats_HM_clusters_km6_pediatric_only_V3.Rdata")
annotation_all = annotation_all[which(annotation_all$Cancer %in% Cancer),]
filtered.norm.RNAseqData = filtered.norm.RNAseqData[,which(colnames(filtered.norm.RNAseqData) %in% annotation_all$Sample)]

load(("./tools/ICR_genes.RData"))
genes_to_correlate = ICR_genes

# Create folders

dir.create(paste0("./Figures/after_split/Correlation_plots/", selected_genes, "_Correlation_plots"), showWarnings = FALSE)

mean_correlation_table = data.frame(Cancertype = Cancer, Mean.correlation = 0)


# Subset RNAseq data to genes to correlate

unavailable.genes <- genes_to_correlate[-which(genes_to_correlate %in% rownames(filtered.norm.RNAseqData))]
subset_RNAseq = t(filtered.norm.RNAseqData[row.names(filtered.norm.RNAseqData) %in% genes_to_correlate, ])
subset_RNAseq_log2 = log(subset_RNAseq +1, 2)

# Correlation matrix
ICR_cor <- cor (subset_RNAseq_log2,method=test)

# Correlation significance
cor.mtest <- function(mat, conf.level = 0.95) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat <- lowCI.mat <- uppCI.mat <- matrix(NA, n, n)
  diag(p.mat) <- 0
  diag(lowCI.mat) <- diag(uppCI.mat) <- 1
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], conf.level = conf.level)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
      lowCI.mat[i, j] <- lowCI.mat[j, i] <- tmp$conf.int[1]
      uppCI.mat[i, j] <- uppCI.mat[j, i] <- tmp$conf.int[2]
    }
  }
  return(list(p.mat, lowCI.mat, uppCI.mat))
}
ICR_cor_sign <- cor.mtest(ICR_cor, 0.95)

# Correlation plot
png(paste0("./Figures/after_split/Correlation_plots/021.", selected_genes, "_Correlation_plots/", selected_genes, "_", test, "_Correlation_plot_", Cancer, ".png"),
    res=600,height=6,width=6,unit="in")


#dev.new()
cex.before <- par("cex")
par(cex = 0.45)
lims=c(-1,1)



if (length(ICR_cor[ICR_cor<0]) == 0) {lims=c(0,1)}
annotation = data.frame (gene = rownames(ICR_cor),color = c(rep("#CC0506",20)),stringsAsFactors = FALSE)
annotation$color[annotation$gene %in% c("IDO1","CD274","CTLA4","FOXP3","PDCD1")] = "#41719C"

annotation = annotation[corrMatOrder(ICR_cor,order="FPC"),]


mean_correlation = round(mean(ICR_cor),2)

corrplot.mixed (ICR_cor,
                #type="lower",
                #p.mat = ICR_cor_sign[[1]],                                                                      # add significance to correlations
                #col = colpattern,
                lower = "square",
                upper ="number",
                order="FPC",
                cl.lim=lims,                                                                                               
                tl.pos ="lt",
                tl.col = as.character(annotation$color),
                insig= "pch",                                                                                              # remove insignificant correlations
                pch = "x",
                pch.cex= 1.5,
                tl.cex = 1/par("cex"),
                cl.cex = 1/par("cex"),
                cex.main = 1/par("cex"),
                mar=c(6,4.1,7,5))
title(main = list(paste0( " Correlation between ", selected_genes, " genes. \n ","Mean: ", mean_correlation,". Number of patients: ", nrow(subset_RNAseq), "."),
                  cex = 2.2), line = -2.5)
title(sub = list(paste0("Figure: TCGAbiolinks normalized, log transformed gene expression data for " ,Cancer, " was \n obtained from TARGET dataset, ",
                        "Significance level of correlation is represented by the size of the squares."), cex = 1.5), line = 1.5)
par(cex = cex.before)
dev.off()


#cat(paste0("For ", Cancer, " mean correlation is ", mean_correlation), file = Log_file, append = TRUE, sep = "\n")
mean_correlation_table$Mean.correlation[mean_correlation_table$Cancertype == Cancer] = mean_correlation


dir.create(paste0("./Analysis/after_split/Correlations"), showWarnings = FALSE)
save(ICR_cor, ICR_cor_sign, file = paste0("./Analysis/after_split/Correlations/Correlation_matrix_ICR_", test, "_", Cancer, ".Rdata"))


#dir.create("./Analysis/Correlations", showWarnings = FALSE)
#save(mean_correlation_table, file = paste0("./Analysis/TCGA_Assembler/Pan_Cancer/Correlation/Correlation_", selected_genes, ".Rdata"))
