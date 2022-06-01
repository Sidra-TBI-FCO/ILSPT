#Script for the event free survival curve per cancer
# Setup environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(master.location)
setwd(paste0(master.location,"/TBI-LAB - Project - Pediatric Pan Cancer TARGET"))
# Install packages and load
source(paste0(toolbox.path,"/R scripts/ipak.function.R"))
required.packages = c("survival","reshape","ggplot2","plyr","Rcpp","colorspace","texreg")
required.bioconductor.packages = "survival"
ipak(required.packages)
ibiopak(required.bioconductor.packages)

#Setting the parameters
ICR_k = "HML_classification"                                                                                            # "HML_classification"
Surv_cutoff_years = 10
Cancer = "OS"

#Loading the required files  
source("~/Sidra Medicine - Research Division/TBI-LAB - General/bioinformatics tools/R scripts/ggkm_Jessica_Pancancer.R")
load(paste0("./Analysis/after_split/ICR_data_OS/TARGET_OS_table_cluster.RData"))
clustering$HML.ICR.Cluster = gsub("ICR High" , "ICR High + Medium",clustering$HML.ICR.Cluster)
clustering$HML.ICR.Cluster = gsub("ICR Medium" , "ICR High + Medium",clustering$HML.ICR.Cluster)

# Create folders
dir.create("./Analysis/after_split/",showWarnings = FALSE)                                                                        # Create folder to save processed data (by Assembler module B)
dir.create(paste0("./Analysis/after_split/Survival_Analysis"), showWarnings = FALSE)

#row.names(clustering) = substring(rownames(clustering),1,16) 
Survival_data =read.csv(paste0("./Data/Clinical_Data/",Cancer,"_Clinical_data.event.free.surv.csv"),stringsAsFactors = FALSE)

rownames(clustering) = gsub("\\." , "-",rownames(clustering))
Survival_data$ICR_cluster = clustering$HML.ICR.Cluster[match(Survival_data$TARGET.USI,rownames(clustering))]
Survival_data$ICRscore = clustering$ICRscore[match(Survival_data$TARGET.USI,rownames(clustering))]

dir.create(paste0("./Figures/after_split/Kaplan_Meier_Plots/aftersplit/",Cancer,"_Event_Free_Survival_KM") , showWarnings = FALSE)

Survival_data = Survival_data[!is.na(Survival_data$ICR_cluster),]
#Survival_data = Survival_data[!duplicated(Survival_data$submitter_id),]


#Survival_data$First.Event[which(Survival_data$Fisrt.Event %in%
#                                               c("Death", "Relapse", "Event", "Progression", "Second Malignant Neoplasm"))] = "Event"
Survival_data$First.Event = gsub("Relapse","Event",Survival_data$First.Event)
Survival_data$First.Event = gsub("Progression","Event",Survival_data$First.Event)
Survival_data$First.Event = gsub("Second Malignant Neoplasm","Event",Survival_data$First.Event)
Survival_data$First.Event = gsub("Death Without Remission","Event",Survival_data$First.Event)
Survival_data$First.Event = gsub("Death","Event",Survival_data$First.Event)


Survival_data$First.Event = gsub("Censored","Event Free",Survival_data$First.Event)
Survival_data$First.Event = gsub("None","Event Free",Survival_data$First.Event)


Y = Surv_cutoff_years * 365
TS.No.Event = Survival_data[Survival_data$First.Event == "Event Free", c("First.Event", "Event.Free.Survival.Time.in.Days", "ICR_cluster","ICRscore")]
colnames(TS.No.Event) = c("First.Event","Time", "Group","ICRscore")
TS.No.Event$Time = as.numeric(as.character(TS.No.Event$Time))
TS.No.Event$Time[TS.No.Event$Time > Y] = Y

#In Case of NBL
#TS.Event = Survival_data[Survival_data$First.Event == "Event", c("First.Event", "Event.Free.Survival.Time.in.Days", "ICR_cluster")]
TS.Event = Survival_data[Survival_data$First.Event == "Event", c("First.Event", "Event.Free.Survival.Time.in.Days", "ICR_cluster","ICRscore")]

colnames(TS.Event) = c("First.Event","Time", "Group","ICRscore")
TS.Event$Time = as.numeric(as.character(TS.Event$Time))
TS.Event$First.Event[which(TS.Event$Time> Y)] = "Event Free"
TS.Event$Time[TS.Event$Time > Y] = Y

TS.Surv = rbind (TS.Event,TS.No.Event)
TS.Surv$Time = as.numeric(as.character(TS.Surv$Time))
TS.Surv$First.Event <- TS.Surv$First.Event == "Event"
TS.Surv = subset(TS.Surv,TS.Surv$Time > 1)                                                                                         # remove patients with less then 1 day follow up time

TS.Surv[,"Group"] = factor(TS.Surv[,"Group"], levels = c("ICR High + Medium", "ICR Low"))

# survival curve
msurv = Surv(TS.Surv$Time/30.4, TS.Surv$First.Event)                                                                                    # calculate the number of months
mfit = survfit(msurv~TS.Surv$Group,conf.type = "log-log")

# Calculations
mdiff = survdiff(eval(mfit$call$formula), data = eval(mfit$call$data))
pval = pchisq(mdiff$chisq,length(mdiff$n) - 1,lower.tail = FALSE)
pvaltxt = ifelse(pval < 0.0001,"p < 0.0001",paste("p =", signif(pval, 3)))


#TS.Surv[,"Group"] = as.factor(TS.Surv[,"Group"])
mHR = coxph(formula = msurv ~ TS.Surv[,"Group"],data = TS.Surv, subset = TS.Surv$Group %in% c("ICR High + Medium", "ICR Low"))
mHR.extract = extract(mHR, include.aic = TRUE,
                            include.rsquared = TRUE, include.maxrs=TRUE,
                            include.events = TRUE, include.nobs = TRUE,
                            include.missings = TRUE, include.zph = TRUE)
HRtxt = paste("Hazard-ratio =", signif(exp(mHR.extract@coef),3),"for",names(mHR$coefficients))
beta = coef(mHR)
se   = sqrt(diag(mHR$var))
p    = 1 - pchisq((beta/se)^2, 1)
CI   = confint(mHR)
CI   = round(exp(CI),2)
PLOT_P = round(p[1],5)
PLOT_HR = round(signif(exp(mHR.extract@coef),2)[1], 2)
PLOT_CI1 = CI[1,1]

PLOT_CI2 = CI[1,2]


dir.create(paste0("./Figures/after_split/Kaplan_Meier_Plots/",Cancer,"_Event_Free_Survival_KM"), showWarnings = FALSE)

# plots
png(paste0("./Figures/after_split/Kaplan_Meier_Plots/",Cancer,"_Event_Free_Survival_KM/",Cancer,"_H_L_Continous_Kaplan_Meier_at_cutoff",Surv_cutoff_years,"_years.png"),
    res=600,height=6,width=8,unit="in")                                                                                           # set filename
dev.new()
ggkm(mfit,
     timeby=12,
     ystratalabs = levels(TS.Surv[,"Group"]),
     ystrataname = NULL,
     main= paste0("Event Free Survival curve across ICR High+ Medium and Low in" , Cancer, "Tumor at cutoff" ,Surv_cutoff_years),
     xlabs = "Time in months",
     palette = c("orange", "blue"),
     PLOT_HR = PLOT_HR,
     PLOT_P = PLOT_P,
     PLOT_CI1 = PLOT_CI1,
     PLOT_CI2 = PLOT_CI2)
dev.off()

# Cox regression (continuous)
uni_variate_ICRscore = coxph(formula = Surv(Time, First.Event) ~ ICRscore, data = TS.Surv)
summary(uni_variate_ICRscore)
