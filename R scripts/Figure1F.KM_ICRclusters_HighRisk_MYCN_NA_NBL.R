# Script for calculating and plotting the Overall survival using the ggsurvplot 

#Setup environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(master.location)
setwd(paste0(master.location,"/TBI-LAB - Project - Pediatric Pan Cancer TARGET"))
# Install packages and load
source(paste0(toolbox.path,"/R scripts/ipak.function.R"))
required.packages = c("survival","survminer")
ipak(required.packages)

## Set Parameters
group.of.interest = "risk_group"
Cancer = "NBL"
Surv_cutoff_years = 10

# Load data
load("./Analysis/after_split/Signature_Enrichment/Thorsson_clustering/Neuroblastoma_clinical_clusters.Rdata")
load(paste0("./Analysis/after_split/ICR_data_",Cancer,"/TARGET_",Cancer,"_table_cluster.RData")) #Pediatric only 

dir.create("./Figures/after_split", showWarnings = FALSE)
dir.create(paste0("./Figures/after_split/Boxplots/NBL/Boxplots_ICR_by_",group.of.interest), showWarnings = FALSE)

clinical_NBL = annotation
clinical_NBL = clinical_NBL[!is.na(clinical_NBL$risk_group),]
clinical_NBL<-clinical_NBL[!(clinical_NBL$risk_group)=="Unknown",]
clinical_NBL<-clinical_NBL[!(clinical_NBL$mycn_status)=="Unknown",]

clinical_NBL= clinical_NBL[which(clinical_NBL$risk_group %in% "High Risk"),]
clinical_NBL= clinical_NBL[which(clinical_NBL$mycn_status %in% "Not Amplified"),]

load(paste0("./Data/Clinical_Data/all_pediatric_only_patients_clinical_data.Rdata"))
Survival_data = clinical_all

Survival_data$ICR_cluster = clustering$HML.ICR.Cluster[match(Survival_data$submitter_id,rownames(clustering))]
Survival_data$ICRscore = clustering$ICRscore[match(Survival_data$submitter_id,rownames(clustering))]

Survival_data$risk_group = annotation$risk_group[match(Survival_data$submitter_id,annotation$Sample)]
Survival_data$mycn_status = annotation$mycn_status[match(Survival_data$submitter_id,annotation$Sample)]
Survival_data = Survival_data[!is.na(Survival_data$risk_group),]
Survival_data = Survival_data[!is.na(Survival_data$mycn_status),]
Survival_data = Survival_data[!duplicated(Survival_data$submitter_id),]
Survival_data = Survival_data[!is.na(Survival_data$ICR_cluster),]

Survival_data<-Survival_data[which(Survival_data$submitter_id %in% clinical_NBL$Sample),]

dir.create(paste0("./Figures/after_split/Kaplan_Meier_Plots/",Cancer,"_Overallsurv_KM") , showWarnings = FALSE)
Y = Surv_cutoff_years * 365
TS.Alive = Survival_data[Survival_data$vital_status == "Alive", c("vital_status", "days_to_last_follow_up", "ICR_cluster","ICRscore")]
colnames(TS.Alive) = c("Status","Time", "Group","ICRscore")
TS.Alive$Time = as.numeric(as.character(TS.Alive$Time))
TS.Alive$Time[TS.Alive$Time > Y] = Y
TS.Dead = Survival_data[Survival_data$vital_status == "Dead", c("vital_status", "days_to_last_follow_up", "ICR_cluster","ICRscore")]
colnames(TS.Dead) = c("Status","Time", "Group","ICRscore")
TS.Dead$Time = as.numeric(as.character(TS.Dead$Time))
TS.Dead$Status[which(TS.Dead$Time> Y)] = "Alive"
TS.Dead$Time[TS.Dead$Time > Y] = Y

TS.Surv = rbind (TS.Dead,TS.Alive)
TS.Surv$Time = as.numeric(as.character(TS.Surv$Time))
TS.Surv$Status <- TS.Surv$Status == "Dead"
TS.Surv = subset(TS.Surv,TS.Surv$Time > 1)                                                                                         # remove patients with less then 1 day follow up time

TS.Surv[,"Group"] = factor(TS.Surv[,"Group"], levels = c("ICR High","ICR Medium","ICR Low"))

# survival curve
dev.new()
fit = survfit(Surv(TS.Surv$Time/30.4, TS.Surv$Status) ~ TS.Surv[,"Group"], data = TS.Surv)

png(paste0("././Figures/after_split/Kaplan_Meier_Plots/",Cancer,"_Overallsurv/036." , Cancer,"_high_risk_mycn.NA_Kaplan_Meier_at_cutoff_",Surv_cutoff_years,"_years_ggsurvplot.png"),
    res=600,height=6,width=6,unit="in")

sur = ggsurvplot(fit,
                 # legend = c(0.2, 0.2),
                 break.time.by = 12, 
                # palette = c("#DF536B","S2"="#7ECD60","S3"="#4A93DF","S4"="#6DDFE3","S5"="#BC3DB5","S6"="#FFA500"),
                 palette = c("ICR High"="red","ICR Medium"="green","ICR Low" ="blue"), 
                 # conf.int = TRUE, 
                # pval.method = TRUE,
                 legend.labs = c("ICR High","ICR Medium","ICR Low"),
                 #legend = "right",
                 pval = TRUE,
                 risk.table = TRUE, 
                 risk.table.y.text.col = TRUE ,
                 #risk.table.col = "strata",
                 #fun = "cumhaz",
                 ylab = "Overall survival probability",
                 fun = function(y) y*100,
                 font.main = 18,
                 font.x =  16,
                 font.y = 16,
                 font.tickslab = 10,
                #tables.theme = theme_cleantable(),  # Clean theme for tables
              
                 #  risk.table.height = 0.2   #Height between groups 
)
sur$table <- sur$table + theme(axis.line = element_blank())
sur$plot <- sur$plot + labs(title = "Survival Curve")

print(sur)

dev.off()
