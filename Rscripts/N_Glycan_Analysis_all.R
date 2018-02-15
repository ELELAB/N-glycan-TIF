 
# Author: Thilde Bagger Terkelsen
# Contact: thilde@cancer.dk
# Place of employment: Computational Biology Laboratory, Danish Cancer Society Research Center, Copenhagen, Denmark
# Date 29-05-2017




# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------





											### LOADING AND PREPARING DATA ###





# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# LOAD DATA AND METADATA
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Glycan data
require(openxlsx)
TIFNIF <- read.xlsx("TIFNIF_all.xlsx", colNames = TRUE, rowNames = TRUE)

# Metadata
TIFNIFinfo <- read.xlsx("TIFNIF_all_info.xlsx", colNames = TRUE, rowNames = TRUE)



# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# REMOVE OUTLIER SAMPELS, TECHNICAL REPLICATES AND SAMPLES WITH LOW TP (< 40%)
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


# Samples 118 and 123 are outliers and were removed

remove <- c(which(TIFNIFinfo$tp %in% c("10", "30", "40")), which(TIFNIFinfo$ID %in% c("TIF66C", "TIF81a", "TIF81b","TIF109.1", "NIF109.1", "TIF118", "NIF118", "TIF123", "TIF123")), which(TIFNIFinfo$Tumor_subtype_corrected_2015_11_20 == "Apocrine"))

# Remove samples from data and metadata
TIFNIF <- TIFNIF[, -remove]
TIFNIFinfo <- TIFNIFinfo[-remove,]

# log2 transformation
TIFNIF <- log2(TIFNIF)

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------







											### PRELIMINARY ANALYSIS ###









# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# GET VECTORS WITH SUBTYPE, CONDITION, BATCH AND NUMBER OF TUMOR SAMPLES
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Condition tumor and normal
NC <- as.factor(TIFNIFinfo$NC)

# BC subtypes
TS <- as.factor(TIFNIFinfo$Tumor_subtype_corrected_2015_11_20)

# Samples batch
batch <- as.factor(TIFNIFinfo$plate)

# Number of Tumor samples
Tn <- length(grep("cancer", TIFNIFinfo$NC))

# Color vectors for plotting Tumor subtypes (TS.cols) and normal+cancer (NC.cols)
TS.cols <- c("purple", "violetred4", "violetred", "orchid2", "orange", "grey60" )
NC.cols <- c("red2", "grey50")



# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# BATCH_CORRECTION CANCER AND NORMAL - PLOTTING ONLY
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Combat batch correction - not used directly for DAA but later for plotting

# design matrix
mod_design <-  model.matrix(~NC)
# batch correction with ComBat
batch_corr_NC <- ComBat(TIFNIF, batch, mod_design, par.prior=TRUE,prior.plots=FALSE)

# Multidimensional scaling plot before and after batch correction, colored by condition
myMDSplot(TIFNIF, NC, NC, NC.cols)
myMDSplot(batch_corr_NC, NC, NC, NC.cols)




# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# BATCH_CORRECTION SUBTYPES AND NORMAL  - PLOTTING ONLY
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# design matrix
mod_design <-  model.matrix(~TS)
# batch correction with ComBat
batch_corr_TS <- ComBat(TIFNIF, batch, mod_design, par.prior=TRUE,prior.plots=FALSE)

# Multidimensional scaling plot before and after batch correction, colored by BC subtype
myMDSplot(TIFNIF, TS, TS, TS.cols)
myMDSplot(batch_corr_TS, TS, TS, TS.cols)



# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------






											### DIFFERENTIAL ABUNDANCE ANALYSIS ###





# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# DIFFERENTIAL ABUNDANCE ANALYSIS (DAA) - CANCER vs. NORMAL
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Model matrix
mod_design <-  model.matrix(~0+NC)

# Contrast for LIMMA
cotr_N_C <- makeContrasts("NCcancer-NCnormal", levels=mod_design)

# DAA, calling function DA_glycans, logFC cutoff 0 and FDR < 0.05
glycans_NC <- DA_glycan(cotr_N_C, TIFNIF, mod_design, 0, 0.05)




# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# DIFFERENTIAL ABUNDANCE ANALYSIS - TUMOR SUBTYPES and NORMAL
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Create all combinations of BC subtypes for LIMMA DA analysis
combinations <- data.frame(t(combn(paste0("TS", levels(TS)), 2)))
combinations$contr <- apply(combinations[,colnames(combinations)], 1, paste, collapse = "-")


# Design incorporating batch
mod_design <-  model.matrix(~0+TS+batch)

# Making group contrasts 
contrast.matrix <- eval(as.call(c(as.symbol("makeContrasts"),as.list(as.character(combinations$contr)),levels=list(mod_design))))


# DAA all comparisons, calling function DA_glycan_apply

# No block. Subtype + batch in design matrix
glycans_TS_1 <- DA_glycan_apply(contrast.matrix, TIFNIF, mod_design, 0, 0.05, NULL, FALSE)

# No block. Subtype in design matrix
mod_design <- model.matrix(~0+TS)
glycans_TS_2 <- DA_glycan_apply(contrast.matrix[1:length(levels(TS)),], TIFNIF, mod_design, 0, 0.05, NULL, FALSE)

# Block with batch. Subtype in design matrix
glycans_TS_3 <- DA_glycan_apply(contrast.matrix[1:length(levels(TS)),], TIFNIF, mod_design, 0, 0.05, batch, FALSE)




# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# DIFFERENTIAL ABUNDANCE ANALYSIS (DAA) - SUBTYPES 
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Create all combinations of BC subtypes for LIMMA DA analysis
combinations <- data.frame(t(combn(paste0("TS_NN", levels(TS_NN)), 2)))
combinations$contr <- apply(combinations[,colnames(combinations)], 1, paste, collapse = "-")


# Design incorporating batch
mod_design <-  model.matrix(~0+TS_NN+batch_NN)

# Making group contrasts 
contrast.matrix <- eval(as.call(c(as.symbol("makeContrasts"),as.list(as.character(combinations$contr)),levels=list(mod_design))))


# DAA all comparisons, calling function DA_glycan_apply
glycans_TS_NN_1 <- DA_glycan_apply(contrast.matrix, TIFNIF[,1:Tn], mod_design, 0, 0.05, NULL, FALSE)

# No block. Subtype in design matrix
mod_design <- model.matrix(~0+TS_NN)
glycans_TS_NN_2 <- DA_glycan_apply(contrast.matrix[1:length(levels(TS_NN)),], TIFNIF[,1:Tn], mod_design, 0, 0.05, NULL, FALSE)

# Block with batch. Subtype in design matrix
glycans_TS_NN_3 <- DA_glycan_apply(contrast.matrix[1:length(levels(TS_NN)),], TIFNIF[,1:Tn], mod_design, 0, 0.05, batch_NN, FALSE)




# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------






							###  TUMOR INFILTRATING LYMPHOCYTE STATUS AND N-GLYCAN ABUNDANCE ###






# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# DATASET - TIF SAMPLES ONLY
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# remove samples with no TIL status 
remove <- unique(c(which(is.na(TIFinfo$TILS)),  which(TIFinfo$TILS %in% c("T3_outside_tumor", "T0_8T0_stroma"))))

# TIL status vector
TILS <- as.factor(as.character(TIFinfo$TILS))

# Simple TIL status, high, low.
TILS <- ifelse(TILS %in% c("T0", "T1"), "low", as.character(TILS))
TILS <- ifelse(TILS %in% c("T2", "T3"), "high", as.character(TILS))

# N-glycans DA between samples with high and low TILs
resTILS <- DA_all_contrasts(TIF, TILS, 0, 0.05, c("high", "low"), remove)



# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# DAA INDIVIDUAL TILS - CD45, CD3, CD4, CD8, CD68
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# remove samples with no CD45 status 
remove <- which(is.na(TIFinfo$CD45))

# CD45 status vector
CD45 <- factor(as.character(TIFinfo$CD45))

resCD45 <- DA_all_contrasts(TIF, CD45, 0, 0.05, c("CD45three", "CD45two", "CD45one"), remove)


# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# remove samples with no CD3 status 
remove <- which(is.na(TIFinfo$CD3))

# CD3 status vector
CD3 <- factor(as.character(TIFinfo$CD3))

resCD3 <- DA_all_contrasts(TIF, CD3, 0, 0.05, c("CD3three", "CD3two", "CD3one"), remove)


# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# remove samples with no CD4 status 
remove <- which(is.na(TIFinfo$CD4))

# CD4 status vector
CD4 <- factor(as.character(TIFinfo$CD4))

resCD4 <- DA_all_contrasts(TIF, CD4, 0, 0.05, c("CD4three", "CD4two", "CD4one"), remove)


# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# remove samples with no CD8 status
remove <- sort(c(which(is.na(TIFinfo$CD8)), which(TIFinfo$CD8 == "CD8three")))

# CD8 status vector
CD8 <- factor(as.character(TIFinfo$CD8))

resCD8 <- DA_all_contrasts(TIF, CD8, 0, 0.05, c("CD8two", "CD8one"), remove)


# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# remove samples with no CD68 status
remove <- which(is.na(TIFinfo$CD68))

# CD68 status vector
CD68 <- factor(as.character(TIFinfo$CD68))

resCD68 <- DA_all_contrasts(TIF, CD68, 0, 0.05, c("CD68three", "CD68two", "CD68one"), remove)



# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------






											### CLINICAL PARAMETERS AND N_GLYCAN ABUNDANCE ###







# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# DAA - Clinical Variables
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Estrogen Receptor

# vector with ER status
ER <- ifelse(TIFinfo$ER == "ER+", "ERp", "ERm")
ER <- as.factor(ER)

resER <- DA_all_contrasts(TIF, ER, 0, 0.05, c("ERp", "ERm"))


# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Progesterone receptor

# vector with PGR status
PGR <- ifelse(TIFinfo$PGR == "PGR+", "PGRp", "PGRm")
PGR <- as.factor(PGR)

resPGR <- DA_all_contrasts(TIF, PGR, 0, 0.05, c("PGRp", "PGRm"))


# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Androgen receptor

# remove samples with no AR status
remove <- which(is.na(TIFinfo$AR))

# vector with AR status
AR <- ifelse(TIFinfo$AR == "AR+", "ARp", "ARm")
AR <- as.factor(AR)

resAR <- DA_all_contrasts(TIF, AR, 0, 0.05, c("ARp", "ARm"), remove)


# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# HER2 receptor

# vector with HER2 status
HER2 <- gsub("[+]", "", TIFinfo$HER2)
HER2 <- as.factor(as.character(paste0(rep("H", length(HER2)), HER2)))

resHER2 <- DA_all_contrasts(TIF, HER2, 0, 0.05, c("H3", "H2", "H1", "H0"))


# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Tumor Grade

# remove samples with no tumor grade
remove <- which(is.na(TIFinfo$Gr))

# vector with tumor grade
GR <- gsub("[+]", "", TIFinfo$Gr)
GR <- as.factor(as.character(paste0(rep("G", length(GR)), GR)))

resGR <- DA_all_contrasts(TIF, GR, 0, 0.05, c("G3", "G2", "G1"), remove)

