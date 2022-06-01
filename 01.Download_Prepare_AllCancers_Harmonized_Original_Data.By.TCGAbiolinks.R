#Setting working environment 
rm(list = ls())
setwd("~/Dropbox (TBI-Lab)/DB-LAB/Projects - Data/02 Active Projects/SDR400151_PIFR_2020_WH_Pediatric_Cancer_TARGET")

#Loading the required packages 
source("./tools/ipak.function.R")
required.packages = c("TCGAbiolinks")
ibiopak(required.packages)
library(TCGAbiolinks)
#Loading the required files
load("~/Dropbox (TBI-Lab)/DB-LAB/Projects - Data/02 Active Projects/SDR400151_PIFR_2020_WH_Pediatric_Cancer_TARGET/tools/geneInfo.July2017.RData")

#Setting parameters

Cancer1 = "TARGET-RT"
Cancer = "RT"

# Creating the query 
query <- GDCquery(project = Cancer1,
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "HTSeq - Counts")

#Downloading the gene expression quantification data and save it 
GDCdownload(query, method = "client", directory = "./Data")

#Preparing the data (creating the Matrix)
# The SummarizedExperiment container contains one or more assays, each represented by a matrix-like object of numeric or other mode. The rows typically represent genomic ranges of interest and the columns represent samples.
GDCprepare(query, save =TRUE ,
           save.filename = paste0("./Processed_Data/001_Exp_",Cancer,"_Processed_Data.rda"), 
           directory = "./Data",
           summarizedExperiment = FALSE)

#Preparing the genes ensembel Id in the data file by removing the numbers after 15 , (removing the ensembel version number)
#length(unique(data$X1))
data$X1 <- substr(data$X1,1,15)
data$X2 <- geneInfo$hgnc_symbol[match(data$X1,geneInfo$ensembl_gene_id)]
data <- data[!is.na(data$X2),]
rownames(data) <- data$X2

data$X1 <- NULL
data$X2 <- NULL
data <- as.matrix(data)
mode(data)
dim(data)
duplicated(colnames(data))


#In case of (RT) To remove the outlier sample 
#test = data.frame(sample = colnames(data), mean = colMeans(data))
#data= as.data.frame(data)
#data$`TARGET-52-PASXGF-01A-01R`= NULL


save(data,file = paste0("./Processed_Data/001_Exp_",Cancer,"_Processed_Data.rda"))

#Create a boxplot correlation and AAIC plot to define outliers 
#WTRnaseq_CorOutliers <- TCGAanalyze_Preprocessing(data)
                                              


