#Script for the Overall survival table (HR table) for all Cancers

# Setup environment
rm(list = ls())
setwd("~/Sidra Medicine - Research Division/TBI-LAB - Project - Pediatric Pan Cancer TARGET/")

# Install packages and load
source("~/Sidra Medicine - Research Division/TBI-LAB - Project - Pediatric Pan Cancer TARGET/tools/ipak.function.R")
required.packages = c("survival","reshape","ggplot2","plyr","Rcpp","colorspace","texreg")
required.bioconductor.packages = "survival"
ipak(required.packages)
ibiopak(required.bioconductor.packages)

#Setting the parameters
ICR_k = "ICRscore"                                                                                            # "HML_classification"
Surv_cutoff_years = 10
#Cancer = "WT"

#Loading the required files  
source("~/Sidra Medicine - Research Division/TBI-LAB - General/bioinformatics tools/R scripts/ggkm_Jessica_Pancancer.R")
#load("./tools/ICR_genes.RData")

TARGET_datasets = read.csv("~/Sidra Medicine - Research Division/TBI-LAB - Project - Pediatric Pan Cancer TARGET/Data/TARGET.dataset.csv", stringsAsFactors = FALSE)
All_Cancers = TARGET_datasets$cancerType

HR_table = data.frame(Cancer = All_Cancers, p_value = 0, HR = 0, CI_lower = 0, CI_upper = 0)

i=6
for (i in 1:length(All_Cancers)){
  Cancer = All_Cancers[i]
  load(paste0("~/Sidra Medicine - Research Division/TBI-LAB - Project - Pediatric Pan Cancer TARGET/Analysis/after_split/ICR_data_",Cancer,"/TARGET_",Cancer,"_table_cluster.RData"))
  
  # Create folders
  dir.create("./Analysis",showWarnings = FALSE)                                                                 
  dir.create(paste0("~/Sidra Medicine - Research Division/TBI-LAB - Project - Pediatric Pan Cancer TARGET/Analysis/Survival_Analysis"), showWarnings = FALSE)
  
  Survival_data=read.csv(paste0("~/Sidra Medicine - Research Division/TBI-LAB - Project - Pediatric Pan Cancer TARGET/Data/Clinical_Data/TARGET-",Cancer,"_clinical.csv"))
  Survival_data$submitter_id = as.character(Survival_data$submitter_id)
  Survival_data  = Survival_data[!is.na(Survival_data$days_to_last_follow_up),]
  Survival_data$ICRscore = clustering$ICRscore[match(Survival_data$submitter_id,rownames(clustering))]
  
  dir.create(paste0("~/Sidra Medicine - Research Division/TBI-LAB - Project - Pediatric Pan Cancer TARGET/Figures/after_split/Kaplan_Meier_Plots/",Cancer,"_Overallsurv_KM"), showWarnings = FALSE)
  
  Survival_data = Survival_data[!is.na(Survival_data$ICRscore),]
  
  Y = Surv_cutoff_years * 365
  TS.Alive = Survival_data[Survival_data$vital_status == "Alive", c("vital_status", "days_to_last_follow_up", "ICRscore")]
  colnames(TS.Alive) = c("Status","Time", "Group")
  TS.Alive$Time = as.numeric(as.character(TS.Alive$Time))
  TS.Alive$Time[TS.Alive$Time > Y] = Y
  
  TS.Dead = Survival_data[Survival_data$vital_status == "Dead", c("vital_status", "days_to_last_follow_up", "ICRscore")]
  
  colnames(TS.Dead) = c("Status","Time", "Group")
  TS.Dead$Time = as.numeric(as.character(TS.Dead$Time))
  TS.Dead$Status[which(TS.Dead$Time> Y)] = "Alive"
  TS.Dead$Time[TS.Dead$Time > Y] = Y
  
  TS.Surv = rbind (TS.Dead,TS.Alive)
  TS.Surv$Time = as.numeric(as.character(TS.Surv$Time))
  TS.Surv$Status <- TS.Surv$Status == "Dead"
  TS.Surv = subset(TS.Surv,TS.Surv$Time > 1)                            #remove patients with less then 1 day follow up time
  
  uni_variate = coxph(formula = Surv(Time, Status) ~ get("Group"), data = TS.Surv)
  summary = summary(uni_variate)
  HR = summary$conf.int[1]
  CI_lower = summary$conf.int[3]
  CI_upper = summary$conf.int[4]
  p_value = summary$coefficients[5]
  HR_table$p_value[which(HR_table$Cancer == Cancer)] = p_value
  HR_table$CI_lower[which(HR_table$Cancer == Cancer)] = CI_lower
  HR_table$CI_upper[which(HR_table$Cancer == Cancer)] = CI_upper
  HR_table$HR[which(HR_table$Cancer == Cancer)] = HR
}

save(HR_table, file = "~/Sidra Medicine - Research Division/TBI-LAB - Project - Pediatric Pan Cancer TARGET/Analysis/after_split/Survival_Analysis/040.HR_table_ICR_Score_Per_Cancer.Rdata")
