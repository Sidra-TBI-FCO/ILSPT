# Setup microenviroment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(master.location)
setwd(paste0(master.location,"/TBI-LAB - Project - Pediatric Pan Cancer TARGET"))
# Install packages and load
source(paste0(toolbox.path,"/R scripts/ipak.function.R"))
required.packages = c("stringr", "ggplot2", "ggpubr","reshape2","tibble","dplyr")
ipak(required.packages)

# Set parameters
Cancer = "PanCancer"
test = "kruskal.test"
######## For Boxplots ########

# Load data
CIBERSORT = read.csv("./Analysis/after_split/CIBERSORTx/CIBERSORTx_all_408_patients.csv",stringsAsFactors = FALSE)

CIBERSORT$Mixture = gsub("\\.","-",CIBERSORT$Mixture) 
rownames(CIBERSORT) = CIBERSORT$Mixture

CIBERSORT$Mixture = NULL
CIBERSORT$P.value = NULL
CIBERSORT$Correlation = NULL
CIBERSORT$RMSE = NULL

## To combine the columns into 6 main groups of cells 
CIBERSORT$Lymphocytes = (CIBERSORT$B.cells.naive + CIBERSORT$B.cells.memory + CIBERSORT$Plasma.cells + CIBERSORT$T.cells.CD8 + CIBERSORT$T.cells.CD4.naive+
CIBERSORT$T.cells.CD4.memory.resting+CIBERSORT$T.cells.CD4.memory.activated+CIBERSORT$T.cells.follicular.helper+
CIBERSORT$T.cells.regulatory..Tregs.+CIBERSORT$T.cells.gamma.delta+CIBERSORT$T.cells.gamma.delta +CIBERSORT$NK.cells.resting+
CIBERSORT$NK.cells.activated)
CIBERSORT$Macrophages = (CIBERSORT$Macrophages.M0 + CIBERSORT$Macrophages.M1 +CIBERSORT$Macrophages.M2+CIBERSORT$Monocytes)
CIBERSORT$Denderitc.cells = (CIBERSORT$Dendritic.cells.activated+ CIBERSORT$Dendritic.cells.resting)
CIBERSORT$Mast.cells = (CIBERSORT$Mast.cells.activated+ CIBERSORT$Mast.cells.resting)
CIBERSORT$Neutrophils = (CIBERSORT$Neutrophils)
CIBERSORT$Eosinophils = (CIBERSORT$Eosinophils)
CIBERSORT = CIBERSORT[,which(colnames(CIBERSORT) %in% c("Lymphocytes","Macrophages","Denderitc.cells","Mast.cells","Neutrophils","Eosinophils"))]
##

proportions = t(CIBERSORT)
load("./Analysis/after_split/Signature_Enrichment/Thorsson_clustering/030.Thorsson_samples_cancers_10000repeats_HM_clusters_km6_pediatric_only_V3.Rdata")

# Analysis
data = data.frame(Sample_ID = annotation_all$Sample,
                  clusters = annotation_all$cluster,
                  Pathway_proportions = NA)

data$clusters = factor(data$clusters, levels = c("S1", "S2", "S3","S4","S5","S6"))

dir.create("./Figures/after_split/Boxplots", showWarnings = FALSE)
dir.create(paste0("./Figures/after_split/Boxplots/CiberSortx_sub.types"), showWarnings = FALSE)
dir.create(paste0("./Figures/after_split/Boxplots/CiberSortx_sub.types/Cibersortex_immune_subtypes_BoxPlot"), showWarnings = FALSE)
dir.create(paste0("./Figures/after_split/Boxplots/CiberSortx_sub.types/Cibersortex_immune_subtypes_BoxPlot/",test),showWarnings = FALSE)

proportions = as.matrix(proportions)

dev.new()
# Plot without p-value
i = 5 
for (i in 1:nrow(proportions)){
  pathway = rownames(proportions)[i]
  data$Pathway_proportions = proportions[pathway,][match(data$Sample_ID, colnames(proportions))]
  plot = ggplot(data, aes(x = clusters, y = Pathway_proportions, fill = clusters)) +
    geom_boxplot(outlier.shape=NA) +
    geom_jitter(size = 0.1, width = 0.15) +
    theme_bw() +
    scale_fill_manual(values = c("#DF536B","#7ECD60","#4A93DF","#6DDFE3","#BC3DB5","#FFA500")) +
    ylab(paste0(pathway)) +
    xlab("") +
   # stat_compare_means(method = "kruskal.test", comparisons = list(c("S1", "S2","S4","S5","S6")),size =2) +
    theme(axis.text.x = element_text(color = "black", size = 10),
          axis.text.y = element_text(color = "black", size = 10),
          axis.title.x = element_text(color = "black", size = 10),
          axis.title.y = element_text(color = "black", size = 16))
  
  png(file = paste0("./Figures/after_split/Boxplots/CiberSortx_sub.types/Cibersortex_immune_subtypes_BoxPlot/",test,"/053.Combined_6_Cibersortex_clusters_",pathway,"_",test,".png"), 
      res = 600, units = "in", width = 5, height = 4)
  plot(plot)
  dev.off()
}

### Boxplot of Kruskal Walis test 
dev.new()
i = 21 
for (i in 1:nrow(proportions)){
  pathway = rownames(proportions)[i]
  data$Pathway_proportions = proportions[pathway,][match(data$Sample_ID, colnames(proportions))]
  plot = ggplot(data, aes(x = clusters, y = Pathway_proportions, fill = clusters)) +
    geom_boxplot(outlier.shape=NA) +
    geom_jitter(size = 0.1, width = 0.15) +
    theme_bw() +
    scale_fill_manual(values = c("#DF536B","#7ECD60","#4A93DF","#6DDFE3","#BC3DB5","#FFA500")) +
    ylab(paste0(pathway)) +
    xlab("") +
    stat_compare_means(method = "kruskal.test",size =4) +
    theme(axis.text.x = element_text(color = "black", size = 10),
          axis.text.y = element_text(color = "black", size = 10),
          axis.title.x = element_text(color = "black", size = 10),
          axis.title.y = element_text(color = "black", size = 18))
  
  png(file = paste0("./Figures/after_split/Boxplots/CiberSortx_sub.types/Cibersortex_immune_subtypes_BoxPlot/",test,"/053.Cibersortex_clusters_",pathway,"_",test,".png"), 
      res = 600, units = "in", width = 5, height = 4)
  plot(plot)
  dev.off()
}

######## For Barplots ########
rm(list = ls())
load("~/R.Config.Rdata")
setwd(master.location)
setwd(paste0(master.location,"/TBI-LAB - Project - Pediatric Pan Cancer TARGET"))
# Install packages and load
source(paste0(toolbox.path,"/R scripts/ipak.function.R"))
required.packages = c("stringr", "ggplot2", "ggpubr","reshape2","tibble","dplyr")
ipak(required.packages)

# Set parameters
Cancer = "PanCancer"
CIBERSORT = read.csv("./Analysis/after_split/CIBERSORTx/CIBERSORTx_all_408_patients.csv",stringsAsFactors = FALSE)

load("./Analysis/after_split/Signature_Enrichment/Thorsson_clustering/030.Thorsson_samples_cancers_10000repeats_HM_clusters_km6_pediatric_only_V3.Rdata")

rownames(CIBERSORT) = CIBERSORT$Mixture

rownames(CIBERSORT) = gsub("\\.","-",rownames(CIBERSORT))

CIBERSORT$Mixture = NULL
CIBERSORT$P.value = NULL
CIBERSORT$Correlation = NULL
CIBERSORT$RMSE = NULL

## To combine the columns into 6 main groups of cells 
CIBERSORT$Lymphocytes = (CIBERSORT$B.cells.naive + CIBERSORT$B.cells.memory + CIBERSORT$Plasma.cells + CIBERSORT$T.cells.CD8 + CIBERSORT$T.cells.CD4.naive+
                           CIBERSORT$T.cells.CD4.memory.resting+CIBERSORT$T.cells.CD4.memory.activated+CIBERSORT$T.cells.follicular.helper+
                           CIBERSORT$T.cells.regulatory..Tregs.+CIBERSORT$T.cells.gamma.delta+CIBERSORT$T.cells.gamma.delta +CIBERSORT$NK.cells.resting+
                           CIBERSORT$NK.cells.activated)
CIBERSORT$Macrophages = (CIBERSORT$Macrophages.M0 + CIBERSORT$Macrophages.M1 +CIBERSORT$Macrophages.M2+CIBERSORT$Monocytes)
CIBERSORT$Denderitc.cells = (CIBERSORT$Dendritic.cells.activated+ CIBERSORT$Dendritic.cells.resting)
CIBERSORT$Mast.cells = (CIBERSORT$Mast.cells.activated+ CIBERSORT$Mast.cells.resting)
CIBERSORT$Neutrophils = (CIBERSORT$Neutrophils)
CIBERSORT$Eosinophils = (CIBERSORT$Eosinophils)
CIBERSORT = CIBERSORT[,which(colnames(CIBERSORT) %in% c("Lymphocytes","Macrophages","Denderitc.cells","Mast.cells","Neutrophils","Eosinophils"))]
##

CIBERSORT$sampleID =rownames(CIBERSORT)
CIBERSORT$cluster = annotation_all$cluster[match(CIBERSORT$sampleID,annotation_all$Sample)]

CIBERSORT = melt(CIBERSORT,id.vars = c("sampleID","cluster"))

CIBERSORT = CIBERSORT[!is.na(CIBERSORT$cluster),]
# Grouped
color_table <- tibble(
  cluster = c("S1", "S2", "S3","S4","S5","S6"),
  Color = c("#DF536B","#7ECD60","#4A93DF","#6DDFE3","#BC3DB5","#FFA500"))

CIBERSORT$cluster <- factor(CIBERSORT$cluster, levels = color_table$cluster)

dev.new()

png(file = paste0("./Figures/after_split/CIBERSORTx/053.Combined.6.groups.Barchart_Cibersortex_proportion_immune.subtypes.png"), 
    res = 600, units = "in", width = 10, height = 7)

plot =ggplot(CIBERSORT, aes(fill=cluster, y=value, x=variable)) + 
  #geom_jitter(size = 0.1, width = 0.15)+
  ylab(paste0("Proportion")) +
  xlab("CIBERSORTx immune cells") +
  scale_fill_manual(values = color_table$Color,drop = FALSE)+
#  theme(text = element_text(size=20),
 #       axis.text.x = element_text(angle=90, hjust=1))+
  theme(text = element_text(size=20),
        axis.text.x = element_text(angle=90, hjust=1),strip.background = element_blank(),strip.text = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank()
  ) +
  stat_summary(geom = "bar", fun = mean, position = "dodge") +
  stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge")

plot(plot)
dev.off()


#### To calculate the Standard deviation and draw the error bars: 

data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

df2 <- data_summary(CIBERSORT, varname="value", 
                    groupnames=c("variable", "cluster"))
# Convert dose to a factor variable
df2$variable=as.factor(df2$variable)
head(df2)

dev.new()

png(file = paste0("./Figures/after_split/CIBERSORTx/Barchart_Cybersortex_proportion_clusters_pediatric_only_sd_ErrorBars.png"), 
    res = 600, units = "in", width = 18, height = 7)

# Keep only upper error bars
ggplot(df2, aes(x=variable, y=value, fill=cluster)) + 
  geom_bar(stat="identity", color="black", position=position_dodge()) +
  geom_errorbar(aes(ymin=value, ymax=value+sd), width=.2,
                position=position_dodge(.9)) +
  ylab(paste0("Proportion")) +
  xlab("CIBERSORTx immune cells") +
  scale_fill_manual(values = color_table$Color,drop = FALSE)+
  theme(text = element_text(size=20),
        axis.text.x = element_text(angle=90, hjust=1),strip.background = element_blank(),strip.text = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank()
        )
dev.off()
