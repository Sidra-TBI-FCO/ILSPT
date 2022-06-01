#Script for the generating forestplot for immune subtypes pancancer

# Setup environment
rm(list = ls())

load("~/R.Config.Rdata")
setwd(master.location)
setwd(paste0(master.location,"/TBI-LAB - Project - Pediatric Pan Cancer TARGET"))
# Install packages and load
source(paste0(toolbox.path,"/R scripts/ipak.function.R"))
required.packages = c("survival","reshape","ggplot2","plyr","Rcpp","colorspace","texreg","forestplot")

required.bioconductor.packages = "survival"
ipak(required.packages)
ibiopak(required.bioconductor.packages)

#Parameters 
Surv.cutoff.years = 10
# Loading the data

HR_table = read.csv("./Analysis/after_split/Survival_Analysis/HR_Thorsson_for_forest_plot.csv",stringsAsFactors = FALSE)


n_clusters = nrow(HR_table)
x = n_clusters + 2

HR.matrix = as.matrix(HR_table)
rownames(HR.matrix) = HR.matrix[,1]
HR.matrix = HR.matrix[,-c(1)]
mode(HR.matrix) = "numeric"

# Cochrane data from the 'rmeta'-package
cochrane_from_rmeta =
  structure(list(
    mean  = as.numeric(c(NA,HR_table$HR[1:n_clusters]), NA),
    lower = c(NA,HR_table$CI_lower[c(1:n_clusters)], NA),
    upper = c(NA,HR_table$CI_upper[c(1:n_clusters)], NA)),
    Names = c("mean", "lower", "upper"),
    row.names = c(NA, -x),
    class = "data.frame")


HR_table$p_value = signif(HR_table$p_value, 3)

HR_table$HR = signif(HR_table$HR, 3)

tabletext<-cbind(
  c("cluster", as.character(HR_table$X)[c(1:n_clusters)]),
  c("p_value", HR_table$p_value[c(1:n_clusters)]),
  c("HR",      HR_table$HR[c(1:n_clusters)]))

dir.create("./Figures/after_split/Forest_plots", showWarnings = FALSE)
pdf(file = paste0("./Figures/after_split/Forest_plots/Thorss_clusters_Forest/024.ImmuneSubtypes_manually_plus_NBL_subgroups_cutoff_", Surv.cutoff.years,".pdf"),
    height = 5, width = 5)

dev.new()
forestplot(mean = HR.matrix[,"HR"],
           lower = HR.matrix[,"CI_lower"],
           upper = HR.matrix[,"CI_upper"],
           labeltext = tabletext[-1,],
           new_page = FALSE,
           zero = 1,
           #is.summary=c(TRUE,rep(FALSE,n.cells),TRUE,rep(FALSE,n.cells),TRUE,FALSE),
           #    clip=c(0.001,2000),
           xlog=TRUE,
           xlim = c(0, 4 , 8 ,12),
           boxsize = .115,
           vertices = FALSE,
           col=fpColors(box="darkblue",line="darkgrey"),
           txt_gp = fpTxtGp(label = gpar(fontsize = 13), xlab = gpar(fontsize = 13),
                            ticks = gpar(fontsize = 12))
)
dev.off()
