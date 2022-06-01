rm(list = ls())
setwd("~/Dropbox (TBI-Lab)/DB-LAB/Projects - Data/02 Active Projects/SDR400151_PIFR_2020_WH_Pediatric_Cancer_TARGET")
load("~/Dropbox (TBI-Lab)/DB-LAB/Projects - Data/02 Active Projects/SDR400151_PIFR_2020_WH_Pediatric_Cancer_TARGET/tools/geneInfo.July2017.RData")

#Loading the required packages 
source("./tools/ipak.function.R")
required.packages = c("TCGAbiolinks","limma")
ibiopak(required.packages)
source("./tools/TCGA-Assembler-2-master/TCGA-Assembler/Module_A.R")
source("./tools/TCGA-Assembler-2-master/TCGA-Assembler/Module_B.R")
TCGASampleTypeFile = paste0("./tools/TCGA-Assembler-2-master/TCGA-Assembler/SupportingFiles/TCGASampleType.txt")
#Setting the parameters
Cancer = "WT"
#Loading the data 
#load("./Processed_Data/001_Exp_NBL_161_Processed_Data.Rdata")
load(paste0("./Processed_Data/001_Exp_",Cancer,"_Processed_Data.rda"))


#Filtering the Primary Solid tumor (TP)
RNAseqData.to.normalize = ExtractTissueSpecificSamples(inputData = data,
                                                       tissueType = c("TP"),
                                                       singleSampleFlag = FALSE,
                                                       sampleTypeFile = TCGASampleTypeFile)

colnames(RNAseqData.to.normalize)<- substr(colnames(RNAseqData.to.normalize),1,16)
duplicated(colnames(RNAseqData.to.normalize))
RNAseqData.to.normalize= RNAseqData.to.normalize[,!duplicated(colnames(RNAseqData.to.normalize))]
RNAseqData.to.normalize <- as.matrix(RNAseqData.to.normalize)
mode(RNAseqData.to.normalize)
dim(RNAseqData.to.normalize)

save(RNAseqData.to.normalize,geneInfo,
     file= paste0("./Processed_Data/002_",Cancer,"_Data_To_normalize_TP_extracted.Rdata"))
load(paste0("./Processed_Data/002_",Cancer,"_Data_To_normalize_TP_extracted.Rdata"))

#Using limma oackage to the density plot before Normalization 
# should plot log2  

RNAseqData.to.normalize <- as.data.frame(RNAseqData.to.normalize)
plotDensities(log(RNAseqData.to.normalize[,1:118],2))
dev.new()


#Data Normalization using TCGAbiolinks
dataNorm <- TCGAanalyze_Normalization(RNAseqData.to.normalize, geneInfo , method = "gcContent" )



plotDensities(log(dataNorm[,1:118],2))
#Open the plot in a different screen 

filtered.norm.RNAseqData = dataNorm
save(filtered.norm.RNAseqData,geneInfo,
     file= paste0("./Processed_Data/003_",Cancer,"_normalized_TCGAbiolinks_TP_filtered.Rdata"))

#sum(rownames(RNAseqData.to.normalize) %in% geneInfo$hgnc_symbol)

