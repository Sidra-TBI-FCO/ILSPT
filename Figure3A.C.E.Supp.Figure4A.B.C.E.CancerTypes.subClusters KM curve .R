# script for calculating and plotting the Overall survival using the ggsurvplot

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
Surv_cutoff_years = 10
Cancer ="NBL_mycn_Amp"      # NBL_mycn_Namp_High   # NBL_mycn_Amp   #NBL_mycn_Namp_Intermed.Low
cluster = c("S1","S2","S3","S5","S6")

dir.create("./Analysis/after_split/",showWarnings = FALSE)                                                                        # Create folder to save processed data (by Assembler module B)
dir.create(paste0("./Analysis/after_split/Survival_Analysis"), showWarnings = FALSE)
load( "./Analysis/after_split/Signature_Enrichment/Thorsson_clustering/030.Thorsson_samples_cancers_10000repeats_HM_clusters_km6_pediatric_only_V3.Rdata")
annotation = annotation_all
annotation = annotation[which(annotation$Cancer %in% Cancer),]
annotation = annotation[which(annotation$cluster %in% cluster),]

#row.names(clustering) = substring(rownames(clustering),1,16) 
Survival_data=read.csv(paste0("./Data/Clinical_Data/046.TARGET-AllCancers_clinical.plus.NBL.csv"))

Survival_data$cluster = annotation$cluster[match(Survival_data$submitter_id,annotation$Sample)]
Survival_data$cancer = annotation$Cancer[match(Survival_data$submitter_id,annotation$Sample)]

Survival_data = Survival_data[!is.na(Survival_data$cluster),]
Survival_data = Survival_data[!duplicated(Survival_data$submitter_id),]

Y = Surv_cutoff_years * 365
TS.Alive = Survival_data[Survival_data$vital_status == "Alive", c("vital_status", "days_to_last_follow_up", "cluster")]
colnames(TS.Alive) = c("Status","Time","cluster")
TS.Alive$Time = as.numeric(as.character(TS.Alive$Time))
TS.Alive$Time[TS.Alive$Time > Y] = Y

TS.Dead = Survival_data[Survival_data$vital_status == "Dead", c("vital_status", "days_to_last_follow_up", "cluster")]

colnames(TS.Dead) = c("Status","Time","cluster")
TS.Dead$Time = as.numeric(as.character(TS.Dead$Time))

TS.Dead$Status[which(TS.Dead$Time> Y)] = "Alive"
TS.Dead$Time[TS.Dead$Time > Y] = Y

TS.Surv = rbind (TS.Dead,TS.Alive)
TS.Surv$Time = as.numeric(as.character(TS.Surv$Time))
TS.Surv$Status <- TS.Surv$Status == "Dead"
TS.Surv = subset(TS.Surv,TS.Surv$Time > 1)                                                                                         # remove patients with less then 1 day follow up time

TS.Surv[,"cluster"] = factor(TS.Surv[,"cluster"], levels = c("S6","S1","S2","S3","S5"))

# survival curve
#dev.new()
fit = survfit(Surv(TS.Surv$Time/30.4, TS.Surv$Status) ~ TS.Surv[,"cluster"], data = TS.Surv)

univartiate = coxph(formula = Surv(Time, Status) ~ cluster , data = TS.Surv)
summary(univartiate)

univartiate = coxph(formula = Surv(Time, Status) ~ cluster , data = TS.Surv)

res.cox.sum <- summary(univartiate)$coefficients
res.cox.sum = as.data.frame(res.cox.sum)
res.cox.sum$HR = round(exp(coef(univartiate)), 2)
res.cox.sum$CI = round(exp(confint(univartiate)), 2)
write.csv(res.cox.sum, file = paste0("./Analysis/after_split/Survival_Analysis/",Cancer,"_for_forest_plot_manually.csv"),quote = FALSE)

#png(paste0("././Figures/after_split/Kaplan_Meier_Plots/",Cancer,"_Overallsurv/" , Cancer,"_thorss_5_Kaplan_Meier_at_cutoff_",Surv_cutoff_years,"_years_HM__km5_10000repeats_pediatric_only_V2_less_colors_ggsurvplot.png"),
#   res=600,height=6,width=6,unit="in")
dir.create(paste0("./Figures/after_split/Kaplan_Meier_Plots/",Cancer,"_Overallsurv"),showWarnings = FALSE)
png(paste0("././Figures/after_split/Kaplan_Meier_Plots/",Cancer,"_Overallsurv/021.",Cancer,"_thorss_5_Kaplan_Meier_at_cutoff_",Surv_cutoff_years,"_years_V3.png"),
    res=600,height=7,width=7,unit="in")   


dev.new()
sur = ggsurvplot(fit,
                 # legend = c(0.2, 0.2),
                 break.time.by = 12, 
                 palette = c("S1"="#DF536B","S2"="#7ECD60","S3"="#4A93DF","S4"="#6DDFE3","S5"="#BC3DB5","S6"="#FFA500"),
                 # palette = c("ICR High"="red","ICR Low" ="blue"), 
                 # conf.int = TRUE, 
                 # pval.method = TRUE,
                 legend.labs = c("S1","S2","S3","S5","S6"),
                 #legend = "left",
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
