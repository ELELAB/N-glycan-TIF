# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Copyright (C) 2018, Thilde Bagger Terkelsen <thilde@cancer.dk> <thildebate@gmail.com>
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------




											### LOADING AND PREPARING DATA ###





# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# LOAD DATA AND METADATA
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Glycan data
TIFNIF <- openxlsx::read.xlsx("TIFNIF_all.xlsx", colNames = TRUE, rowNames = TRUE)

# Metadata
TIFNIFinfo <- openxlsx::read.xlsx("TIFNIF_all_info.xlsx", colNames = TRUE, rowNames = TRUE)


# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# REMOVE SAMPLES WITH TP LOWER THAN 30%
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

remove <- which(TIFNIFinfo$tp %in% c("10"))

# Remove samples from data and metadata
TIFNIF <- TIFNIF[, -remove]
TIFNIFinfo <- TIFNIFinfo[-remove,]


# Only TIF samples
TIF <- TIFNIF[,1:length(grep("cancer", TIFNIFinfo$NC))]
TIFinfo <- TIFNIFinfo[1:length(grep("cancer", TIFNIFinfo$NC)),]


# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# LOAD DATA AND SURVIVAL DATA
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


# Survival data
survivaldata <- read.table("survivaldata.txt", header = TRUE)


# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# REMOVE RECURRANCES, SAMPLES LACKING INFO - column "Relaps_status" is kept for relaps free survival
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# removal of samples from survival data: TIF123 (second and third entry, reccurance), TIF80 (second entry, recurrance) and TIF113 (error in survival info)
remove_surv <- sort(c(which(survivaldata$ID %in% c("TIF80", "TIF123") & survivaldata$number_of_reccurances %in% c(2,3)), which(survivaldata$ID == "TIF113")))

# Remove samples from data
survivaldata <- survivaldata[-remove_surv,]


# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# MATCHING SURVIVAL DATA AND N-GLYCAN ABUNDANCES
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Get overlap in sample ID
ov <- intersect(as.character(survivaldata$ID), as.character(colnames(TIF)))

# Curate survivaldata
survivaldata <- survivaldata[as.character(survivaldata$ID) %in% ov,]

# Calculate outcome in years
survivaldata$time_to_Outcome_years <- survivaldata$time_to_Outcome_months/12

# Curate N-Glycan data
TIF <- TIF[colnames(TIF) %in% ov]
TIFinfo <- TIFinfo[TIFinfo$ID %in% ov,]

# Tranformation log2
TIFlog2 <- log2(TIF)

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------






											### COX REGRESSION NO CONFOUNDERS ###
									



# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# creating survival object

surv_object <- data.frame(t(TIFlog2), survivaldata$time_to_Outcome_months, survivaldata$Age_at_surgery, survivaldata$Outcome_status, survivaldata$Relaps_status, survivaldata$DFBC)
colnames(surv_object) <- c(rownames(TIFlog2), "time", "age", "outcome", "relaps", "DFBC")

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Check for proprotional hazards
pha <- apply(surv_object[,1:63], 2, function(GP) cox.zph(coxph(Surv(time, outcome) ~ pspline(age, df=2) + GP, data = surv_object))$table[,3])
pha <- colnames(pha)[apply(pha, 2, function(x) any(x<0.05))]
print(paste0("The following features may be in violation of the proportional hazard assumption: ", pha))



# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# OVERALL SURVIVAL
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

surv_age_outcome <- apply(surv_object[,1:63], 2, function(GP) coxph(Surv(time, outcome) ~ pspline(age, df=2) + GP, data = surv_object))
my_survival(surv_age_outcome, "Survival_corrected_Age")

# Reform and save data for Figure 3A
#Figure3Adata <- surv_age_outcome
#save(Figure3Adata, file="Figure3Adata.Rdata")



# Survival curves for overall survival for GP5, GP10, GP23, GP24, GP38 and CoreF (GP62). Reform and save data for Figure 6B.
GPs <- my_survival_curve(TIFlog2, survivaldata, c(5,10,23,24,38,62))
#Figure3Bdata_1 <- list(GPs, TIFlog2)
#save(Figure3Bdata_1, file="Figure3Bdata_1.Rdata")


fitSurv5 <- survfit(GPs[[1]], newdata=data.frame(age = c(66, 66), GP = c(mean(c(quantile(as.numeric(TIFlog2[5,]))[4], quantile(as.numeric(TIFlog2[5,]))[5])),mean(c(quantile(as.numeric(TIFlog2[5,]))[1], quantile(as.numeric(TIFlog2[5,]))[2])))))
fitSurv10 <- survfit(GPs[[2]], newdata=data.frame(age = c(66, 66), GP = c(mean(c(quantile(as.numeric(TIFlog2[10,]))[4], quantile(as.numeric(TIFlog2[10,]))[5])), mean(c(quantile(as.numeric(TIFlog2[10,]))[1], quantile(as.numeric(TIFlog2[10,]))[2])))))
fitSurv23 <- survfit(GPs[[3]], newdata=data.frame(age = c(66, 66), GP = c(mean(c(quantile(as.numeric(TIFlog2[23,]))[4], quantile(as.numeric(TIFlog2[23,]))[5])),mean(c(quantile(as.numeric(TIFlog2[23,]))[1], quantile(as.numeric(TIFlog2[23,]))[2])))))
fitSurv24 <- survfit(GPs[[4]], newdata=data.frame(age = c(66, 66), GP = c(mean(c(quantile(as.numeric(TIFlog2[24,]))[4], quantile(as.numeric(TIFlog2[24,]))[5])),mean(c(quantile(as.numeric(TIFlog2[24,]))[1], quantile(as.numeric(TIFlog2[24,]))[2])))))
fitSurv38 <- survfit(GPs[[5]], newdata=data.frame(age = c(66, 66), GP = c(mean(c(quantile(as.numeric(TIFlog2[38,]))[4], quantile(as.numeric(TIFlog2[38,]))[5])),mean(c(quantile(as.numeric(TIFlog2[38,]))[1], quantile(as.numeric(TIFlog2[38,]))[2])))))
fitSurv62 <- survfit(GPs[[6]], newdata=data.frame(age = c(66, 66), GP = c(mean(c(quantile(as.numeric(TIFlog2[62,]))[4], quantile(as.numeric(TIFlog2[62,]))[5])),mean(c(quantile(as.numeric(TIFlog2[62,]))[1], quantile(as.numeric(TIFlog2[62,]))[2])))))
#Figure3Bdata_2 <- list(fitSurv5, fitSurv10, fitSurv23, fitSurv24, fitSurv38, fitSurv62)
#save(Figure3Bdata_2, file="Figure3Bdata_2.Rdata")



# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# SURVIVAL WITH CONFOUNDERS
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Samples with no TIL status, unique status or no info on grade.
remove <- c(which(is.na(TIFinfo$TILS)), which(TIFinfo$TILS == "T3_outside_tumor"), which(is.na(TIFinfo$Gr)))

# Remove samples from datasets
TILSinfo_GR <- TIFinfo[-remove,]
TILS_GR_surv_object <- surv_object[-remove,]

			  
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# OVERALL SURVIVAL WITH CONFOUNDERS
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
			  
# TILS as confounder
surv_age_TILS_outcome <- apply(TILS_GR_surv_object[,1:63], 2, function(GP) coxph(Surv(time, outcome) ~ pspline(age, df=2) + as.factor(TILSinfo_GR$TILS) + GP, data = TILS_GR_surv_object))
my_survival(surv_age_TILS_outcome, "Survival_corrected_Age_TILS")

# Tumor grade as confounder
surv_age_GR_outcome <- apply(TILS_GR_surv_object[,1:63], 2, function(GP) coxph(Surv(time, outcome) ~ pspline(age, df=2) + as.factor(TILSinfo_GR$Gr) + GP, data = TILS_GR_surv_object))
my_survival(surv_age_GR_outcome, "Survival_corrected_Age_Grade")

# TILs + Tumor grade as confounders
surv_age_TILS_GR_outcome <- apply(TILS_GR_surv_object[,1:63], 2, function(GP) coxph(Surv(time, outcome) ~ pspline(age, df=2) + as.factor(TILSinfo_GR$TILS) + as.factor(TILSinfo_GR$Gr) + GP, data = TILS_GR_surv_object))
my_survival(surv_age_TILS_GR_outcome, "Survival_corrected_Age_TILS_Grade")





# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# RELAPS-FREE SURVIVAL
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


surv_age_relaps <- apply(surv_object[,1:63], 2, function(GP) coxph(Surv(time, relaps) ~ pspline(age, df=2) + GP, data = surv_object))
my_survival(surv_age_relaps, "Relapsfree_survival_Age")



# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# RELAPS-FREE SURVIVAL WITH CONFOUNDERS
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


# TILS as confounder
surv_age_TILS_relaps <- apply(TILS_GR_surv_object[,1:63], 2, function(GP) coxph(Surv(time, relaps) ~ pspline(age, df=2) + as.factor(TILSinfo_GR$TILS) + GP, data = TILS_GR_surv_object))
my_survival(surv_age_TILS_relaps, "Relapsfree_corrected_Age_TILS")

# Tumor grade as confounder
surv_age_GR_relaps<- apply(TILS_GR_surv_object[,1:63], 2, function(GP) coxph(Surv(time, relaps) ~ pspline(age, df=2) + as.factor(TILSinfo_GR$Gr) + GP, data = TILS_GR_surv_object))
my_survival(surv_age_GR_relaps, "Relapsfree_corrected_Age_Grade")

# TILs + Tumor grade as confounders
surv_age_TILS_GR_relaps <- apply(TILS_GR_surv_object[,1:63], 2, function(GP) coxph(Surv(time, relaps) ~ pspline(age, df=2) + as.factor(TILSinfo_GR$TILS) + as.factor(TILSinfo_GR$Gr) + GP, data = TILS_GR_surv_object))
my_survival(surv_age_TILS_GR_relaps, "Relapsfree_corrected_Age_TILS_Grade")

				 
				 

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# CAUSE OF DEATH KNOW TO BE BREAST CANCER
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


surv_age_DFBC <- apply(surv_object[,1:63], 2, function(GP) coxph(Surv(time, DFBC) ~ pspline(age, df=2) + GP, data = surv_object))
my_survival(surv_age_DFBC, "DFBCsurvival_Age")

#Supfig2data <- surv_age_DFBC
#save(Supfig2data, file="Supfig2data.Rdata")

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# CAUSE OF DEATH KNOW TO BE BREAST CANCER WITH CONFOUNDERS
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


# TILS as confounder
surv_age_TILS_DFBC <- apply(TILS_GR_surv_object[,1:63], 2, function(GP) coxph(Surv(time, DFBC) ~ pspline(age, df=2) + as.factor(TILSinfo_GR$TILS) + GP, data = TILS_GR_surv_object))
my_survival(surv_age_TILS_DFBC, "DFBCsurvival_corrected_Age_TILS")

# Tumor grade as confounder
surv_age_GR_DFBC<- apply(TILS_GR_surv_object[,1:63], 2, function(GP) coxph(Surv(time, DFBC) ~ pspline(age, df=2) + as.factor(TILSinfo_GR$Gr) + GP, data = TILS_GR_surv_object))
my_survival(surv_age_GR_DFBC, "DFBCsurvival_corrected_Age_Grade")

# TILs + Tumor grade as confounders
surv_age_TILS_GR_DFBC <- apply(TILS_GR_surv_object[,1:63], 2, function(GP) coxph(Surv(time, DFBC) ~ pspline(age, df=2) + as.factor(TILSinfo_GR$TILS) + as.factor(TILSinfo_GR$Gr) + GP, data = TILS_GR_surv_object))
my_survival(surv_age_TILS_GR_DFBC, "DFBCsurvival_corrected_Age_TILS_Grade")
