# Setup environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(master.location)
setwd(paste0(master.location,"/TBI-LAB - Project - Pediatric Pan Cancer TARGET"))

# Install packages and load
source(paste0(toolbox.path,"/R scripts/ipak.function.R"))
required.packages <- c("corrplot", "stringr")
ipak(required.packages)                                                                                                  # Install and load required packages

# Set Parameters
Cancer = "PanCancer"                                                                                                    # Specify download method (this information to be used when saving the file)
Cancer_skip = c("")
colpattern = colorRampPalette(c("#152B7E", "white", "#1B7E09"))(n = 430)
Gene.set = "Thorsson"
test = "spearman"

#Loading the data

load(paste0("./Analysis/after_split/Signature_Enrichment/GSEA_", Cancer, 
            "_pediatric_only", Gene.set,".Rdata"))


# Create folders
dir.create("./Figures/after_split/",showWarnings = FALSE)
dir.create(paste0("./Figures/after_split/Correlation_plots"), showWarnings = FALSE)
dir.create(paste0("./Figures/after_split/Correlation_plots/", Gene.set, "_Correlation_plots"), showWarnings = FALSE)

Hallmark.enrichment.score.df = as.data.frame(t(ES))
#Hallmark.enrichment.score.df$ICR_score = clustering$ICRscore[match(row.names(Hallmark.enrichment.score.df), row.names(clustering))]

Hallmark_GSEA_cor <- cor (Hallmark.enrichment.score.df,method=test)

mean_correlation_table = data.frame(Cancertype = Cancer, Mean.correlation = 0)


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
Hallmark_GSEA_cor_sign <- cor.mtest(Hallmark_GSEA_cor, 0.95)

# Correlation plot
png(paste0("./Figures/after_split/Correlation_plots/", Gene.set, "_Correlation_plots/028", Gene.set, "_", test, "_Correlation_plot_", Cancer, "_pediatric_only_hclust_test.png"),
    res=400,height=20,width=20,unit="in")
#dev.new()
cex.before <- par("cex")
par(cex = 0.6)
lims=c(-1,1)
if (length(Hallmark_GSEA_cor[Hallmark_GSEA_cor<0]) == 0) {lims=c(0,1)}

#rownames(Hallmark_GSEA_cor) = gsub(".*] ", "",rownames(Hallmark_GSEA_cor))


annotation = data.frame (gene = rownames(Hallmark_GSEA_cor),color = c(rep("#CC0506",105)),stringsAsFactors = FALSE)
#annotation$color[annotation$gene]
annotation = annotation[corrMatOrder(Hallmark_GSEA_cor,order="hclust"),]
mean_correlation = round(mean(Hallmark_GSEA_cor),2)

dev.new()
corrplot.mixed (Hallmark_GSEA_cor,
                #type="lower",
                #p.mat = ICR_cor_sign[[1]],                                                                      # add significance to correlations
                #col = colpattern,
                addrect = 5,
                lower = "square",
                upper ="square",
                order="hclust",
                cl.lim=lims,                                                                                               # only positive correlations
                tl.pos ="n",    #n Or "lt" for the titles 
                tl.col = as.character(annotation$color),
                insig= "pch",                                                                                              # remove insignificant correlations
                pch = "x",
                pch.cex= 4,
                tl.cex = 0.8/par("cex"),
                cl.cex = 1/par("cex"),
                cex.main = 1/par("cex"),
                mar=c(6,4.1,7,5))

title(main = list(paste0( " Correlation between ", Gene.set, " genes. \n ","Mean: ", mean_correlation,"."),
                  cex = 4), line = -6)
#title(sub = list(paste0("Figure: TCGAbiolinks normalized, log transformed gene expression data for" ,Cancer, " was \n obtained from TARGET, using ", download.method ,
#                       "Significance level of correlation is represented by the size of the squares."), cex = 1.5), line = 1.5)
par(cex = cex.before)
dev.off()

mean_correlation_table$Mean.correlation[mean_correlation_table$Cancertype == Cancer] = mean_correlation

dir.create(paste0("./Analysis/after_split/Correlations"), showWarnings = FALSE)


save(Hallmark_GSEA_cor, Hallmark_GSEA_cor_sign, file = paste0("./Analysis/after_split/Correlations/029.Correlation_matrix_",Gene.set,"_", test, "_", Cancer, "_pediatric_only.Rdata"))

write.csv(Hallmark_GSEA_cor,file = "./Analysis/after_split/Signature_Enrichment/Thorsson_clustering/029.Correlation_Thorsson_105_signatures_TARGET_paper_supp.csv")
########Create a heatmap for the correlation values 
library(ComplexHeatmap)

Heatmap(Hallmark_GSEA_cor, clustering_distance_rows = "euclidean",clustering_distance_columns = "euclidean",cluster_rows = TRUE ,cluster_columns = TRUE , row_names_max_width = unit(1, "in"),
        show_column_names = TRUE,show_row_names = TRUE, name = "Enrichment\n z score"
)

