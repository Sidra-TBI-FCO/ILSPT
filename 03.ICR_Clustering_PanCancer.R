#Script for the ICR clustering of TARGET PanCancer 

# Setup environment
rm(list = ls())
setwd("~/Sidra Medicine - Research Division/TBI-LAB - Project - Pediatric Pan Cancer TARGET")

# Install packages and load
source("./tools/ipak.function.R")
required.packages = c("ConsensusClusterPlus","clue","heatmap3")
ipak(required.packages)

# Set parameters
groups = "HML.ICR.Cluster" # "k4" or "HL_k4" or "k2" or "HL_k2"
Cancer = "Pancancer"
# load data
load("./Processed_Data/003_PanCancer_pediatric_only_normalized_TCGAbiolinks_TP_filtered.Rdata")
load("./Analysis/after_split/annotation_with_NBL_subgroups.Rdata")
#annotation_all= annotation_all[which(annotation_all$Cancer %in% Cancer),]
filtered.norm.RNAseqData = filtered.norm.RNAseqData[,which(colnames(filtered.norm.RNAseqData) %in% annotation_all$Sample)]

load("~/Sidra Medicine - Research Division/TBI-LAB - General/bioinformatics tools/ICR genes/ICR_genes.RData")

# Create directories
dir.create("./Analysis", showWarnings = FALSE)
dir.create(paste0("./Analysis/after_split/ICR_data_",Cancer), showWarnings = FALSE)
# subset expression matix for ICR
ICR_subset_RNAseq = t(filtered.norm.RNAseqData[which(row.names(filtered.norm.RNAseqData) %in% ICR_genes), ])  #t=transpose :columns will be the rows and rows will be the columns # which rownames are in the ICR genes (the position)
# taking all the columns(the samples), and choosing the rows that have the ICR genes and then transpose.  
ICR_subset_RNAseq_log2 = log(ICR_subset_RNAseq +1, 2)

# Hierarchical Clustering
source("./tools/stefanofunctions.R")
setwd(paste0("./Analysis/after_split/ICR_data_",Cancer))

ddist = dist(ICR_subset_RNAseq_log2)  #making the cluster dendrogram (distance between samples)
class(ddist)
ConsensusClusterObject <- ConsensusClusterPlus(ddist,                                                                
                                               maxK = 6,                                                                              # set K
                                               pItem = 0.8,
                                               reps=5000,                                                                             # set repeats
                                               title=paste0("renormalized.Full.matrix.ICR.reps5000"),              # Output filename (no path)
                                               clusterAlg = "hc",                                                                     # clustering Algorithm : Hierarchiocal clustering
                                               innerLinkage = "ward.D2",                                                              # for color coding the clusters use tmyPal = ...
                                               finalLinkage = "complete",
                                               plot = 'pdf',                                                                          # write resut to pdf (Alt:png)
                                               writeTable = TRUE,
                                               verbose = TRUE)
outputfiles = list.files(paste0("renormalized.Full.matrix.ICR.reps5000"), full.names = TRUE)
class_files = outputfiles[grep("consensusClass", outputfiles)]  #grep :Finds the files that have the (consensusClass) in the name

N.files = length(class_files)  #length: counts the number of elements in the character vector (class_files)
table_cluster_assignment = data.frame(ICRscore = rowMeans(ICR_subset_RNAseq_log2))  # the mean expressions of the rows(ICR genes) is the ICR score, create a dataframe.

j=1
for (j in 1:N.files){
  file = paste0("./", class_files[j])
  consensus_class = read.csv(file = file,header=FALSE)     #load the correct file with the information about the cluster assignment 
  group = paste0("Group_k",j+1)
  colnames(consensus_class) = c("PatientID", group)   #change the column name to patients ID and group
  rownames(consensus_class) = consensus_class$PatientID  # change the row name to patient ID 
  consensus_class$PatientID = NULL      #delete patient ID column
  table_cluster_assignment[,group] = consensus_class[,group][match(rownames(table_cluster_assignment), rownames(consensus_class))] #[,group] : to create a column called group.
  
  transl_table_ICR_cluster = aggregate(ICRscore~get(group),data = table_cluster_assignment, FUN=mean) # calculate the mean ICR score per group
  colnames(transl_table_ICR_cluster) = c(group,"mean_ICRscore")
  transl_table_ICR_cluster = cbind(transl_table_ICR_cluster[order(transl_table_ICR_cluster$mean_ICRscore),],ICR_name=paste0("ICR",c(1:(j+1))))
  
  ICR_cluster = paste0("ICR_cluster_k",j+1)
  table_cluster_assignment[,ICR_cluster] = transl_table_ICR_cluster$ICR_name[match(table_cluster_assignment[,group],
                                                                                   transl_table_ICR_cluster[,group])]
}

#calinsky    # to see how clean is your clusters
sHc <- hclust(ddist, method = "ward.D2")
aCalinsky <- calinsky(sHc,gMax=10)
pdf(file = paste0("./renormalized.Full.matrix.ICR.reps5000/ICR_cluster_assignment_k2-6.Calinsky.pdf"), width = 16, height = 6)
plot(aCalinsky, type = "l", col = "grey", main = "Calinsky & Harabasz curve", xlab = "# of groups")
text(1:length(aCalinsky), aCalinsky, paste(1:length(aCalinsky)))
dev.off()
optimal.calinsky = which(aCalinsky == max(aCalinsky[3:5]))


#save data
save(table_cluster_assignment,optimal.calinsky, file = paste0("./TARGET_",Cancer,"_table_cluster.Rdata")) 

#setwd("../..")

# HML assignment
table_cluster_assignment$HML.ICR.Cluster = NA
table_cluster_assignment$HML.ICR.Cluster[which(table_cluster_assignment$ICR_cluster_k3 == "ICR1")] = "ICR Low"
table_cluster_assignment$HML.ICR.Cluster[which(table_cluster_assignment$ICR_cluster_k3 == "ICR2")] = "ICR Medium"
table_cluster_assignment$HML.ICR.Cluster[which(table_cluster_assignment$ICR_cluster_k3 == "ICR3")] = "ICR High"

clustering = table_cluster_assignment

save(clustering,optimal.calinsky, file =paste0("./TARGET_",Cancer,"_table_cluster.Rdata"))

######################################################################## END OF SCRIPT ################################################

##ICR heatmap3
dir.create("./Figures/Heatmaps", showWarnings = FALSE)
dir.create("./Figures/Heatmaps/ICR_heatmaps", showWarnings = FALSE)

load("./Analysis/ICR_data_NBL/TARGET_NBL_table_cluster.RData")
clustering <- clustering[order(clustering$ICRscore),]
expression.matrix.subset.ICR <- filtered.norm.RNAseqData[ICR_genes,rownames(clustering)]
colors = clustering[,c("ICRscore", "HML.ICR.Cluster")]

if(groups == "HML.ICR.Cluster"){
  colors$HML.ICR.Cluster[colors$HML.ICR.Cluster=="ICR Low"]<-"blue"
  colors$HML.ICR.Cluster[colors$HML.ICR.Cluster=="ICR Medium"]<-"green"
  colors$HML.ICR.Cluster[colors$HML.ICR.Cluster=="ICR High"]<-"red"
}

colors <- colors[,2,drop=FALSE]
colors = as.matrix(colors)

# reoder expression matrix based on colors
ICR_subset_RNAseq_log2 = t(ICR_subset_RNAseq_log2)
expression.matrix.subset.ICR = ICR_subset_RNAseq_log2[rev(ICR_genes),rownames(colors)]

png(filename = paste0("./Figures/Heatmaps/ICR_heatmaps/ICR_heatmapv1_", groups, ".png"), res = 600,
    width = 4, height = 4.5, units = "in")
heatmap3(expression.matrix.subset.ICR,
         main="TCGA COAD ICR RNASeq",
         ColSideColors=colors,
         Colv= NA , 
         Rowv = NA,
         #as.dendrogram(sHc),
         #col=bluered(75) ,
         #labRow = NA,
         #labCol = samples,
         scale='row',
         margins = c(12, 7),
         labCol = NA)

dev.off()
