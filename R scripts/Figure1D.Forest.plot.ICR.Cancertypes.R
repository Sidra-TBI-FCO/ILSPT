# Script for generating the forestplot of ICRscores

# Setup environment
rm(list = ls())

# Install packages and load
required.packages <- c("RColorBrewer", "forestplot")
library(required.packages)

# Set parameters

Surv.cutoff.years = 10
Survival_outcome = "OS"
#Cancer = "PanCancer"
# Load data
load(paste0("./Analysis/after_split/Survival_Analysis/040.HR_table_ICR_Score_Per_Cancer.Rdata"))

HR_table = HR_table[order(HR_table$HR),]

## Forest plot seperate script
n_cancers = nrow(HR_table)
x = n_cancers + 2

HR.matrix = as.matrix(HR_table)
rownames(HR.matrix) = HR.matrix[,1]
HR.matrix = HR.matrix[,-c(1)]
mode(HR.matrix) = "numeric"

#HR_table = HR_table[order(HR_table$HR),]

# Cochrane data from the 'rmeta'-package
cochrane_from_rmeta =
  structure(list(
    mean  = as.numeric(c(NA,HR_table$HR[1:n_cancers]), NA),
    lower = c(NA,HR_table$CI_lower[c(1:n_cancers)], NA),
    upper = c(NA,HR_table$CI_upper[c(1:n_cancers)], NA)),
    Names = c("mean", "lower", "upper"),
    row.names = c(NA, -x),
    class = "data.frame")


HR_table$p_value = signif(HR_table$p_value, 3)
HR_table$HR = signif(HR_table$HR, 3)
tabletext<-cbind(
  c("Cancer", as.character(HR_table$Cancer)[c(1:n_cancers)]),
  c("p-value", HR_table$p_value[c(1:n_cancers)]),
  c("HR",      HR_table$HR[c(1:n_cancers)]))

pdf(file = paste0("./Figures/after_split/Forest_plots/Forest_plot_PanCancer_cutoff_", Surv.cutoff.years,".pdf"),
    height = 7, width = 8)

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
           boxsize = .1,
           vertices = FALSE,
           col=fpColors(box="darkblue",line="darkgrey"),
           txt_gp = fpTxtGp(label = gpar(fontsize = 13), xlab = gpar(fontsize = 13),
                            ticks = gpar(fontsize = 12))
)
dev.off()
