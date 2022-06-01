# Four scripts for generating figures 4A,B,C

######### 1 Forest plot  #############
#Setup environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(master.location)
setwd(paste0(master.location,"/TBI-LAB - Project - Pediatric Pan Cancer TARGET"))
# Install packages and load
source(paste0(toolbox.path,"/R scripts/ipak.function.R"))
required.packages <- c("stringr", "survival","RColorBrewer", "forestplot")
ipak(required.packages)

# Set parameters
Cancer = "PanCancer"
Gene.set = "checkpoints"
Surv.cutoff.years = 10
#cluster= "S4"

# Load data
clinical_data = read.csv(paste0("./Data/Clinical_Data/046.TARGET-AllCancers_clinical.plus.NBL.csv"))
#load(paste0("./Analysis/after_split/Signature_Enrichment/GSEA_",Cancer,"_",Gene.set,".Rdata"))

load("./Processed_Data/003_PanCancer_pediatric_only_normalized_TCGAbiolinks_TP_filtered.Rdata")

load("./Analysis/after_split/Signature_Enrichment/Thorsson_clustering/030.Thorsson_samples_cancers_10000repeats_HM_clusters_km6_pediatric_only_V3.Rdata")
#annotation = annotation[which(annotation$cluster %in% cluster),]
filtered.norm.RNAseqData = filtered.norm.RNAseqData[,which(colnames(filtered.norm.RNAseqData) %in% annotation_all$Sample)]
clinical_data = clinical_data[which(clinical_data$submitter_id %in% annotation_all$Sample),]

#Checkpoints = read.csv("~/Sidra Medicine - Research Division/Sidra Medicine - Research Division/TBI-LAB - General/bioinformatics tools/Checkpoint genes/checkpoint.genes.csv",stringsAsFactors = FALSE)
load("./Analysis/after_split/Signature_Enrichment/checkpoints_1_genes_order.Rdata")
filtered.norm.RNAseqData = filtered.norm.RNAseqData[which(rownames(filtered.norm.RNAseqData) %in% order_signatures),]
ICR_subset_RNAseq_log2 = log(filtered.norm.RNAseqData +1, 2)

clinical_data = clinical_data [!duplicated(clinical_data$submitter_id),]

rownames(clinical_data) = clinical_data$submitter_id

ES = ICR_subset_RNAseq_log2

ES = t(ES)
clinical_data = merge(clinical_data, ES, by = "row.names")

for (i in 1:ncol(ES)){
  col = colnames(ES)[i]
  ES[, col] = (ES[, col] - min(ES[, col]))/(max(ES[,col])-min(ES[,col]))
}

HR_table = data.frame(Signature = colnames(ES), p_value = NA, HR = NA, CI_lower = NA, CI_upper = NA)

i=1
for (i in 1:ncol(ES)){
  Group.of.interest = colnames(ES)[i]
  Y = Surv.cutoff.years * 365
  # time / event object creation
  TS.Alive = clinical_data[clinical_data$vital_status == "Alive", c("vital_status", "days_to_last_follow_up", Group.of.interest)]
  colnames(TS.Alive) = c("Status","Time", Group.of.interest)
  TS.Alive$Time = as.numeric(as.character(TS.Alive$Time))
  TS.Alive$Time[TS.Alive$Time > Y] = Y
  
  TS.Dead = clinical_data[clinical_data$vital_status == "Dead", c("vital_status", "days_to_last_follow_up", Group.of.interest)]
  colnames(TS.Dead) = c("Status","Time", Group.of.interest)
  TS.Dead$Time = as.numeric(as.character(TS.Dead$Time))
  TS.Dead$Status[which(TS.Dead$Time> Y)] = "Alive"
  TS.Dead$Time[TS.Dead$Time > Y] = Y
  
  TS.Surv = rbind (TS.Dead,TS.Alive)
  TS.Surv$Time = as.numeric(as.character(TS.Surv$Time))
  TS.Surv$Status <- TS.Surv$Status == "Dead"
  TS.Surv = subset(TS.Surv,TS.Surv$Time > 1)                               # remove patients with less then 1 day follow up time
  
  uni_variate = coxph(formula = Surv(Time, Status) ~ get(Group.of.interest), data = TS.Surv)
  summary = summary(uni_variate)
  HR = summary$conf.int[1]
  CI_lower = summary$conf.int[3]
  CI_upper = summary$conf.int[4]
  p_value = summary$coefficients[5]
  HR_table$p_value[which(HR_table$Signature == Group.of.interest)] = p_value
  HR_table$CI_lower[which(HR_table$Signature == Group.of.interest)] = CI_lower
  HR_table$CI_upper[which(HR_table$Signature == Group.of.interest)] = CI_upper
  HR_table$HR[which(HR_table$Signature == Group.of.interest)] = HR
}

write.csv(HR_table,file = paste0("./Analysis/after_split/Survival_Analysis/Figure.4A.HR_table_",Cancer,"_ES_",Gene.set,"_cutoff_", Surv.cutoff.years,".csv"))

dir.create("./Analysis/after_split/Survival_Analysis", showWarnings = FALSE)
#save(HR_table, file = paste0("./Analysis/after_split/Survival_Analysis/HR_table_",Cancer,"_ES_",Gene.set,"_cutoff_", Surv.cutoff.years,".Rdata"))

#Create the forest plot
#load(paste0("./Analysis/after_split/Survival_Analysis/HR_table_",Cancer,"_ES_",Gene.set,"_cutoff_", Surv.cutoff.years,".Rdata"))
HR_table = HR_table[order(HR_table$HR),]


## Forest plot seperate script
n_Signature = nrow(HR_table)

x = n_Signature + 2

HR.matrix = as.matrix(HR_table)
rownames(HR.matrix) = HR.matrix[,1]
HR.matrix = HR.matrix[,-c(1)]
mode(HR.matrix) = "numeric"

HR_table = HR_table[order(HR_table$HR),]

# Cochrane data from the 'rmeta'-package
cochrane_from_rmeta =
  structure(list(
    mean  = as.numeric(c(NA,HR_table$HR[1:n_Signature]), NA),
    lower = c(NA,HR_table$CI_lower[c(1:n_Signature)], NA),
    upper = c(NA,HR_table$CI_upper[c(1:n_Signature)], NA)),
    Names = c("mean", "lower", "upper"),
    row.names = c(NA, -x),
    class = "data.frame")


HR_table$p_value = signif(HR_table$p_value, 3)
HR_table$HR = signif(HR_table$HR, 3)
tabletext<-cbind(
  c("Signature", as.character(HR_table$Signature)[c(1:n_Signature)]),
  c("p-value", HR_table$p_value[c(1:n_Signature)]),
  c("HR",      HR_table$HR[c(1:n_Signature)]))

genes_order = rownames(HR.matrix)
save(genes_order,file = "./Analysis/after_split/Signature_Enrichment/checkpoints_2_genes_order.Rdata")

dir.create("./Figures/after_split/Forest_plots/checkpoints_forest.plots", showWarnings = FALSE)
pdf(file = paste0("./Figures/after_split/Forest_plots/checkpoints_forest.plots/checkpoints_across_subtypes ",Cancer,"_Forest_plot_ES",Gene.set,"_cutoff_", Surv.cutoff.years, ".pdf"),
    height = 5, width = 4)
dev.new()
forestplot(mean = HR.matrix[,"HR"],
           lower = HR.matrix[,"CI_lower"],
           upper = HR.matrix[,"CI_upper"],
           labeltext = tabletext[-1,],
           new_page = FALSE,
           zero = 1,
           #is.summary=c(TRUE,rep(FALSE,n.cells),TRUE,rep(FALSE,n.cells),TRUE,FALSE),
           clip=c(0.001,2000),
           xlog=TRUE,
           xlim = c(0, 4 , 8 ,12),
           boxsize = .25,
           vertices = FALSE,
           col=fpColors(box="darkblue",line="darkgrey"),
           txt_gp = fpTxtGp(label = gpar(fontsize = 7), xlab = gpar(fontsize = 7),
                            ticks = gpar(fontsize = 10))
)
dev.off()

######################### 2 Dotted plot  #############################
#DOTTED PLOT OF CLUSTERS 
# Setup environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(master.location)

setwd(paste0(master.location,"/TBI-LAB - Project - Pediatric Pan Cancer TARGET"))

# Install packages and load
source(paste0(toolbox.path,"/R scripts/ipak.function.R"))
required.bioconductor.packages = c("GSVA","ComplexHeatmap", "ggplot2", "ggpubr", "circlize",
                                   "dendsort", "stringr")                                                                   
ibiopak(required.bioconductor.packages)

# Set parameters
Gene.set = "Checkpoints"
Cancer = "PanCancer"

# Load data
#Checkpoints = read.csv("~/Sidra Medicine - Research Division/TBI-LAB - General/bioinformatics tools/Checkpoint genes/checkpoint.genes.csv",stringsAsFactors = FALSE)
checkpoints_inhib_act = read.csv("~/Sidra Medicine - Research Division/TBI-LAB - Project - Pediatric Pan Cancer TARGET/tools/Checkpoints_act_inh.csv",stringsAsFactors = FALSE)
checkpoints_inhib_act = checkpoints_inhib_act[-which(checkpoints_inhib_act$HUGO.GENE.SYM == "KIR3DL3"),]

#Checkpoints$role = checkpoints_inhib_act$Role[match(Checkpoints$HUGO.GENE.SYM,checkpoints_inhib_act$HUGO.GENE.SYM)]

load("./Analysis/after_split/Signature_Enrichment/Thorsson_clustering/030.Thorsson_samples_cancers_10000repeats_HM_clusters_km6_pediatric_only_V3.Rdata")

load(paste0("./Processed_Data/003_PanCancer_pediatric_only_normalized_TCGAbiolinks_TP_filtered.Rdata"))

filtered.norm.RNAseqData = filtered.norm.RNAseqData[which(rownames(filtered.norm.RNAseqData) %in% checkpoints_inhib_act$HUGO.GENE.SYM),]

plot_df = data.frame(t(filtered.norm.RNAseqData))
plot_df$cluster = annotation_all$cluster[match(rownames(plot_df),
                                           annotation_all$Sample)]


#plot_df$Cancer = factor(plot_df$Cancer, levels = c("CCSK", "NBL", "WT" , "OS" , "RT"))
plot_df$cluster = factor(plot_df$cluster, levels = c("S1","S2", "S3", "S4" , "S5" , "S6"))

plot_df_agg = aggregate(.~cluster, plot_df, FUN = median)
rownames(plot_df_agg) = plot_df_agg$cluster
plot_df_agg$cluster = NULL
plot_df_agg = as.matrix(t(plot_df_agg))

# z-score enrichment .matrix
plot_df_agg_z = plot_df_agg
for(j in 1: nrow(plot_df_agg_z))  {
  plot_df_agg_z[j,] = (plot_df_agg[j,]-mean(plot_df_agg[j,]))/sd(plot_df_agg[j,]) # z-score the enrichment matrix
}

plot_df_agg_z = plot_df_agg_z[complete.cases(plot_df_agg_z), ]


load("./Analysis/after_split/Signature_Enrichment/checkpoints_2_genes_order.Rdata")
plot_df_agg_z = plot_df_agg_z[which(rownames(plot_df_agg_z) %in% genes_order),]
checkpoints_inhib_act = checkpoints_inhib_act[which(checkpoints_inhib_act$HUGO.GENE.SYM %in% rownames(plot_df_agg_z)),]
rownames(checkpoints_inhib_act) = checkpoints_inhib_act$HUGO.GENE.SYM
genes_order = genes_order[which(genes_order %in% rownames(plot_df_agg_z))]
checkpoints_inhib_act = checkpoints_inhib_act[genes_order,]
row.names(checkpoints_inhib_act) = NULL


plot_df_agg_z = plot_df_agg_z[checkpoints_inhib_act$HUGO.GENE.SYM,]   #Order the rows to be the same 

write.csv(plot_df_agg_z,file = paste0("./Analysis/after_split/Signature_Enrichment/Figure.4B.zscored.matrix.ES.allCancers.csv"))


ha = HeatmapAnnotation(`cluster` = colnames(plot_df_agg_z),
                       # col = list(`Cancer` = c("CCSK" = "#BB8FCE", "NBL" = "#FFC300", "WT" = "#2980B9" , "OS"= "#1ABC9C" , "RT"="#FA8072")))
                       col = list(`cluster` = c("S1"="#DF536B","S2"="#7ECD60","S3"="#4A93DF","S4"="#6DDFE3","S5"="#BC3DB5","S6"="#FFA500")))

checkpoints_inhib_act$Role = factor(checkpoints_inhib_act$Role, levels = c("activatory", "inhibitory"))


row_ha = rowAnnotation(
  Role = checkpoints_inhib_act$Role,
  col = list(`Role` = c("activatory" = "#D35400", "inhibitory" = "#273746")))

dir.create("./Figures/after_split/Dotted_Heatmap", showWarnings = FALSE)

#save(genes_order,file = "./Analysis/after_split/Signature_Enrichment/checkpoints_genes_order.Rdata")

#dend = dendsort(hclust(dist(plot_df_agg_z)))
dev.new()

#Checkpoints$role = factor(Checkpoints$role, levels = c("activatory","inhibitory"))

#Checkpoints_df <- Checkpoints[ order(Checkpoints$role), ]
#order = Checkpoints_df$HUGO.GENE.SYM
#plot_df_agg_z = plot_df_agg_z[order,]

col_fun = colorRamp2(c(-2,0, 2.16), c("blue","white", "red"))
png(paste0("./Figures/after_split/Dotted_Heatmap/checkoints_clusters_", Gene.set, "_PLUS_ROLE_splitted_pediatric_only.png"),
    res = 600, width = 7, height = 8, units = "in")
dev.new()
Heatmap(plot_df_agg_z, cluster_rows = FALSE, cluster_columns = FALSE,
        show_column_names = FALSE, top_annotation = ha,left_annotation = row_ha ,name = "Mean expression (z-scored)",
        #row_names_gp = gpar(fontsize = 20, fontface = "bold"),
        rect_gp = gpar(type = "none"),
        #row_split = checkpoints_inhib_act$Role ,
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "white", fill = "#EBEBEB"))
          grid.circle(x = x, y = y, r = 0.025,gp = gpar(fill = col_fun(plot_df_agg_z[i, j]), col = NA))
        }
)
dev.off()

##################### 3 calculate the HR tables ######################

#Setup environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(master.location)
setwd(paste0(master.location,"/TBI-LAB - Project - Pediatric Pan Cancer TARGET"))
# Install packages and load
source(paste0(toolbox.path,"/R scripts/ipak.function.R"))
required.packages <- c("stringr", "survival","RColorBrewer", "forestplot",
                       "ComplexHeatmap","ggplot2", "ggpubr","dendsort" ,"circlize")
ipak(required.packages)

# Set parameters
cluster = "S6"
Gene.set = "checkpoints"
#Gene.set2 = "ICR_list"
Surv.cutoff.years = 10
Surv_cutoff_years = 10
#mycn_status = "Not Amplified"
# Load data
#load("./Analysis/after_split/Signature_Enrichment/Thorsson_clustering/Neuroblastoma_clinical_clusters.Rdata")
#annotation_NBL = annotation
#annotation_NBL = annotation_NBL[which(annotation_NBL$mycn_status %in% mycn_status),]

clinical_data = read.csv(paste0("./Data/Clinical_Data/046.TARGET-AllCancers_clinical.plus.NBL.csv"))
clinical_data = clinical_data[!is.na(clinical_data$days_to_last_follow_up),]
load(paste0("./Processed_Data/003_PanCancer_pediatric_only_normalized_TCGAbiolinks_TP_filtered.Rdata"))
checkpoints_inhib_act = read.csv("~/Sidra Medicine - Research Division/TBI-LAB - Project - Pediatric Pan Cancer TARGET/tools/Checkpoints_act_inh.csv",stringsAsFactors = FALSE)
checkpoints_inhib_act = checkpoints_inhib_act[-which(checkpoints_inhib_act$HUGO.GENE.SYM == "KIR3DL3"),]
filtered.norm.RNAseqData = filtered.norm.RNAseqData[which(rownames(filtered.norm.RNAseqData) %in% checkpoints_inhib_act$HUGO.GENE.SYM),]
#load("./Analysis/after_split/Signature_Enrichment/checkpoints_genes_order.Rdata")
#filtered.norm.RNAseqData = filtered.norm.RNAseqData[which(rownames(filtered.norm.RNAseqData) %in% genes_order),]
filtered.norm.RNAseqData = log(filtered.norm.RNAseqData +1, 2)

#save(filtered.norm.RNAseqData,file = "./Analysis/after_split/Signature_Enrichment/checkpoints_pediatric_only_log2transformed_matrix.Rdata")

load("./Analysis/after_split/Signature_Enrichment/Thorsson_clustering/030.Thorsson_samples_cancers_10000repeats_HM_clusters_km6_pediatric_only_V3.Rdata")
annotation = annotation_all
ES = filtered.norm.RNAseqData
annotation = annotation[which(annotation$cluster %in% cluster),]
rownames(clinical_data) = clinical_data$submitter_id
ES = ES[,which(colnames(ES)  %in% annotation$Sample)]

immune_sig_df = as.data.frame(t(ES))

i=1
for (i in 1:ncol(immune_sig_df)){
  col = colnames(immune_sig_df)[i]
  immune_sig_df[, col] = (immune_sig_df[, col] - min(immune_sig_df[, col]))/(max(immune_sig_df[,col])-min(immune_sig_df[,col]))
}

clinical_data = merge(clinical_data, immune_sig_df, by = "row.names")

HR_table = data.frame(Signature = colnames(immune_sig_df), p_value = NA, HR = NA, CI_lower = NA, CI_upper = NA)

i=1
for (i in 1:ncol(immune_sig_df)){
  Group.of.interest = colnames(immune_sig_df)[i]
  Y = Surv.cutoff.years * 365
  # time / event object creation
  TS.Alive = clinical_data[clinical_data$vital_status == "Alive", c("vital_status", "days_to_last_follow_up", Group.of.interest)]
  colnames(TS.Alive) = c("Status","Time", Group.of.interest)
  TS.Alive$Time = as.numeric(as.character(TS.Alive$Time))
  TS.Alive$Time[TS.Alive$Time > Y] = Y
  
  TS.Dead = clinical_data[clinical_data$vital_status == "Dead", c("vital_status", "days_to_last_follow_up", Group.of.interest)]
  colnames(TS.Dead) = c("Status","Time", Group.of.interest)
  TS.Dead$Time = as.numeric(as.character(TS.Dead$Time))
  TS.Dead$Status[which(TS.Dead$Time> Y)] = "Alive"
  TS.Dead$Time[TS.Dead$Time > Y] = Y
  
  TS.Surv = rbind (TS.Dead,TS.Alive)
  TS.Surv$Time = as.numeric(as.character(TS.Surv$Time))
  TS.Surv$Status <- TS.Surv$Status == "Dead"
  TS.Surv = subset(TS.Surv,TS.Surv$Time > 1)                               # remove patients with less then 1 day follow up time
  
  uni_variate = coxph(formula = Surv(Time, Status) ~ get(Group.of.interest), data = TS.Surv)
  summary = summary(uni_variate)
  HR = summary$conf.int[1]
  CI_lower = summary$conf.int[3]
  CI_upper = summary$conf.int[4]
  p_value = summary$coefficients[5]
  HR_table$p_value[which(HR_table$Signature == Group.of.interest)] = p_value
  HR_table$CI_lower[which(HR_table$Signature == Group.of.interest)] = CI_lower
  HR_table$CI_upper[which(HR_table$Signature == Group.of.interest)] = CI_upper
  HR_table$HR[which(HR_table$Signature == Group.of.interest)] = HR
}

dir.create(paste0("./Analysis/after_split/Survival_Analysis/Thorss_clusters_HR/",Gene.set), showWarnings = FALSE)
save(HR_table, file = paste0("./Analysis/after_split/Survival_Analysis/Thorss_clusters_HR/",Gene.set,"/HR_table_",cluster,"_ES_",Gene.set,"_cutoff_", Surv.cutoff.years,"pediatric.Rdata"))


##################### 4 Draw the HR heatmap  ######################

#Setup environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(master.location)
setwd(paste0(master.location,"/TBI-LAB - Project - Pediatric Pan Cancer TARGET"))
# Install packages and load
source(paste0(toolbox.path,"/R scripts/ipak.function.R"))
required.packages <- c("stringr", "survival","RColorBrewer", "forestplot",
                       "ComplexHeatmap","ggplot2", "ggpubr","dendsort" ,"circlize")
ipak(required.packages)

Gene.set = "checkpoints" 
Surv.cutoff.years = 10
Clusters = "ALL"

TARGET.clusters = data.frame(Cluster = c("S1","S2","S3","S4","S5","S6"))

TARGET.clusters$Cluster = as.character(TARGET.clusters$Cluster)

if (Clusters == "ALL") { 
  Clusters = c("Pancancer",TARGET.clusters$Cluster)
}

N.sets = length(Clusters)

load("./Analysis/after_split/Signature_Enrichment/checkpoints_pediatric_only_log2transformed_matrix.Rdata")

ES = filtered.norm.RNAseqData
Signatures = rownames(ES)
df = data.frame(matrix(nrow=34, ncol =7 ))
colnames(df) = Clusters
rownames(df) = Signatures

load("./Analysis/after_split/Signature_Enrichment/checkpoints_2_genes_order.Rdata")
i=1
for (i in 1:N.sets){
  Cluster = Clusters[i]
  load(paste0("./Analysis/after_split/Survival_Analysis/Thorss_clusters_HR/",Gene.set,"/HR_table_",Cluster,"_ES_",Gene.set,"_cutoff_", Surv.cutoff.years,"pediatric.Rdata"))
  HR_table$Signature = as.character(HR_table$Signature)
  rownames(HR_table) = HR_table$Signature
  HR_table = HR_table[genes_order,]
  df[,Cluster] = HR_table$HR[match(rownames(df),HR_table$Signature)]
  df = df[genes_order,]
}

df_HR = data.frame(matrix(nrow=34, ncol = 7))
colnames(df_HR) = Clusters
rownames(df_HR) = Signatures
df_HR = df_HR[genes_order,]

i =1
df_pvalue = df_HR
for (i in 1:N.sets){
  Cluster = Clusters[i]
  load(paste0("./Analysis/after_split/Survival_Analysis/Thorss_clusters_HR/",Gene.set,"/HR_table_",Cluster,"_ES_",Gene.set,"_cutoff_", Surv.cutoff.years,"pediatric.Rdata"))
  rownames(HR_table) = HR_table$Signature
  HR_table = HR_table[genes_order,]
  df_HR[,Cluster] = HR_table$HR[match(rownames(df_HR),HR_table$Signature)]
  df_pvalue[,Cluster] = HR_table$p_value[match(rownames(df_pvalue),HR_table$Signature)]
}

#all_rownames = rownames(df_HR)
#for (i in 1:nrow(df_HR)){
#  rowname = all_rownames[i]
# p_val = df_pvalue[i,]
#if(sum(p_val >= 0.1) == 7 ){
# df_HR = df_HR[-which(rownames(df_HR) == rowname),]
#}
#}

df_HR = as.matrix(df_HR)

df_pvalue_categories = df_pvalue

i=1
j=1
for (i in 1:nrow(df_pvalue_categories)){
  for (j in 1:ncol(df_pvalue_categories)){
    if(df_pvalue_categories[i,j] < 0.05){
      df_pvalue_categories[i,j] = 0.01
    }
    if(df_pvalue_categories[i,j] >= 0.05 & df_pvalue_categories[i,j] < 0.1){
      df_pvalue_categories[i,j] = 0.07
    }
    if(df_pvalue_categories[i,j] >= 0.1){
      df_pvalue_categories[i,j] = 0.3
    }
  }
}

df_pvalue = df_pvalue[rownames(df_HR),]
df_pvalue = -log10(df_pvalue)
df_pvalue = as.matrix(df_pvalue)
df_pvalue_categories =  df_pvalue_categories[rownames(df_pvalue),]

write.csv(df_pvalue,file = paste0("./Analysis/after_split/Survival_Analysis/Figure.4C.pvalue.HR.matrix.allCancers.csv"))
write.csv(df_HR,file = paste0("./Analysis/after_split/Survival_Analysis/Figure.4C.HR.HR.matrix.allCancers.csv"))

ha = HeatmapAnnotation(`Cluster` = colnames(df),
                       col = list(`Cluster` = c("Pancancer" ="brown","S1"="#DF536B","S2"="#7ECD60","S3"="#4A93DF","S4"="#6DDFE3","S5"="#BC3DB5","S6"="#FFA500")
                       ),annotation_name_gp = gpar(fontsize = 22))
dev.new()
col_fun = colorRamp2(c(0,1,2),c("#800000","white","#000080"))
col_fun2 = colorRamp2(c(0.01,0.07,0.3),c("#EBDEF0","#f0dca1","white"))  #Pink color is significant , wellow is not significant


dir.create(paste0("./Figures/after_split/Heatmaps/HR_heatmap/ImmuneSubtypes"),showWarnings = FALSE)
dir.create(paste0("./Figures/after_split/Heatmaps/HR_heatmap/ImmuneSubtypes/",Gene.set),showWarnings = FALSE)
png(filename = paste0("./Figures/after_split/Heatmaps/HR_heatmap/ImmuneSubtypes/",Gene.set,"/","Complexheatmap_", Gene.set, "_PLUS_Pancancer_HR_red_blue_dotted_TEEEEESSSSSST.png"), res = 600,
    width = 8, height = 12, units = "in")


HM_HR = Heatmap(df_HR, cluster_rows = FALSE ,cluster_columns = FALSE , row_names_max_width = unit(6, "in"),column_title_gp = gpar(fontsize = 14),row_names_gp = gpar(fontsize = 22),
                heatmap_legend_param =list(title_gp=gpar(fontsize=10),legend_width=unit(10,"cm")),
                show_column_names = FALSE, top_annotation = ha, name = "HR",column_order = c("Pancancer","S1","S2","S3","S4","S5","S6"),
                rect_gp = gpar(type = "none"),
                cell_fun = function(j, i, x, y, width, height, fill) {
                  grid.rect(x = x, y = y, width = 0.1, height = 0.02, gp = gpar(col = "white", fill = col_fun2(df_pvalue_categories[i,j])))
                  grid.circle(x = x, y = y, r =  abs(df_pvalue[i, j])/2 * min(unit.c(width, height)) ,
                              gp = gpar(fill = col_fun(df_HR[i, j]), col = NA)) 
                })

#  grid.rect(x = x, y = y, width = 0.1, height = 0.015, gp = gpar(col = "white", fill = col_fun2(df_pvalue_categories[i,j])))
draw(HM_HR, heatmap_legend_side = "left", annotation_legend_side = "left")
dev.off()

r.dend <- row_dend(HM_HR)
rcl.list <- row_order(HM_HR)

order_signatures_1 = rownames(df_HR)[rcl.list]

save(order_signatures,file = paste0("./Analysis/after_split/Signature_Enrichment/",Gene.set,"_1_genes_order.Rdata"))
