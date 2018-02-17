# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------




											### LOADING AND PREPARING DATA ###






# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# LOAD DATA AND METADATA
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Glycan data
TIFNIF <-  openxlsx::read.xlsx("paired_TIF_NIF_samples.xlsx", colNames = TRUE, rowNames = TRUE)
#rownames(TIFNIF) <- TIFNIF$X1
#TIFNIF$X1 <- NULL

# Metadata
TIFNIFinfo <-  openxlsx::read.xlsx("sampleinfo_paired_TIF_NIF_samples.xlsx", colNames = TRUE, rowNames = TRUE)

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# REMOVE OUTLIER SAMPELS, TECHNICAL REPLICATES AND SAMPLES WITH LOW TP (< 40%)
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Samples 118 and 123 are outliers and were removed
remove <- c(TIFNIFinfo[which(TIFNIFinfo$tp %in% c("10", "30", "40")),]$patient_ID, TIFNIFinfo[which(TIFNIFinfo$ID %in% c("TIF118", "NIF118", "TIF123", "NIF123")),]$patient_ID, TIFNIFinfo[grep("Apocrine", TIFNIFinfo$Tumor_subtype_corrected_2015_11_20),]$patient_ID)
remove <- which(TIFNIFinfo$patient_ID %in% remove)

# Remove samples from data and metadata
TIFNIF <- TIFNIF[,-remove]
TIFNIFinfo <- TIFNIFinfo[-remove,]

# log2 transformation
TIFNIF <- log2(TIFNIF)

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------







											### PRELIMINARY ANALYSIS ###







# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# GET VECTORS WITH SUBTYPE, CONDITION, BATCH, PATIENT ID AND NUMBER OF TUMOR SAMPLES
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# number of tumor samples
Tn <- length(grep("cancer", TIFNIFinfo$NC))

# Condition tumor and normal
NC <- as.factor(as.character(TIFNIFinfo$NC))

# BC Subtypes
TS <- factor(as.character(TIFNIFinfo$Tumor_subtype_corrected_2015_11_20), levels=c("normal", "LumA", "LumB", "LumB_HER2_enriched", "HER2", "TNBC"))

# BC no Normal
TS_NN <- as.factor(as.character(TS[1:Tn]))

# Patient ID
patient <- as.factor(as.character(TIFNIFinfo$patient))

# Sample batch
batch <- as.factor(TIFNIFinfo$plate)


# Color vectors for plotting Tumor subtypes (TS.cols) and normal+cancer (NC.cols)
TS.cols <- c("grey60", "navyblue", "blue", "darkolivegreen3", "orange", "red")
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
myMDSplot(TIFNIF, NC, "", NC.cols)
myMDSplot(batch_corr_NC, NC, "", NC.cols)




# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# BATCH_CORRECTION SUBTYPES AND NORMAL - PLOTTING ONLY
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# design matrix
mod_design <-  model.matrix(~TS)
# batch correction with ComBat
batch_corr_TS <- ComBat(TIFNIF, batch, mod_design, par.prior=TRUE,prior.plots=FALSE)

# Multidimensional scaling plot before and after batch correction, colored by BC subtype
myMDSplot(TIFNIF, TS, "", TS.cols)
myMDSplot(batch_corr_TS, TS, "", TS.cols)

# Save data for Figure 1
# Figure1Adata <- list(batch_corr_TS, TS, TS.cols)
# save(Figure1Adata, file= "Figure1Adata.Rdata")



# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------






											### DIFFERENTIAL ABUNDANCE ANALYSIS ###





# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# DIFFERENTIAL ABUNDANCE ANALYSIS (DAA) - CANCER vs. NORMAL
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


# Model matrix
mod_design <-  model.matrix(~0+NC+patient)

# Contrast for LIMMA
cotr_N_C <- makeContrasts("NCcancer-NCnormal", levels=mod_design)

# DAA, calling function DA_glycans, logFC cutoff 0 and FDR < 0.05
glycans_NC <- DA_glycan(cotr_N_C, TIFNIF, mod_design, 0, 0.05)

# Reform and save data for Figure 1
#data_ordered <- rbind(glycans_NC[[1]][order(glycans_NC[[1]]$logFC, decreasing=TRUE),], glycans_NC[[2]][order(glycans_NC[[2]]$logFC),])
#Figure1Bdata <- data.frame(rownames(data_ordered), data_ordered$logFC, 1/data_ordered$adj.P.Val, log2(1/data_ordered$adj.P.Val))
#colnames(Figure1Bdata) <- c("GPs", "logFC", "inverseFDR", "scaledinverseFDR")
#save(Figure1Bdata, file="Figure1Bdata.Rdata")  



# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# DIFFERENTIAL ABUNDANCE ANALYSIS (DAA) - SUBTYPES and NORMAL
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


# Create all combinations of BC subtypes for LIMMA DA analysis
combinations <- data.frame(t(combn(paste0("TS", levels(TS)), 2)))
combinations$contr <- apply(combinations[,colnames(combinations)], 1, paste, collapse = "-")


# Design incorporating batch and patient
mod_design <-  model.matrix(~0+TS+patient)

# Making group contrasts 
contrast.matrix <- eval(as.call(c(as.symbol("makeContrasts"),as.list(as.character(combinations$contr)),levels=list(mod_design))))


# DAA all comparisons, calling function DA_glycan_apply

# No block. Subtype + patient in design matrix
glycans_TS_1 <- DA_glycan_apply(contrast.matrix, TIFNIF, mod_design, 0, 0.05, NULL, FALSE)

# Block for batch. Subtype + patient in design matrix
glycans_TS_2 <- DA_glycan_apply(contrast.matrix, TIFNIF, mod_design, 0, 0.05, batch, FALSE)

# Block for patient. Subtype in design matrix
mod_design <- model.matrix(~0+TS)
glycans_TS_3 <- DA_glycan_apply(contrast.matrix[1:length(levels(TS)),], TIFNIF, mod_design, 0, 0.05, patient, FALSE)


# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# DIFFERENTIAL ABUNDANCE ANALYSIS (DAA) - SUBTYPES 
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


# Create all combinations of BC subtypes for LIMMA DA analysis
combinations <- data.frame(t(combn(paste0("TS_NN", levels(TS_NN)), 2)))
combinations$contr <- apply(combinations[,colnames(combinations)], 1, paste, collapse = "-")


# Design incorporating batch and patient
mod_design <-  model.matrix(~0+TS_NN)

# Making group contrasts 
contrast.matrix <- eval(as.call(c(as.symbol("makeContrasts"),as.list(as.character(combinations$contr)),levels=list(mod_design))))


# DAA all comparisons, calling function DA_glycan_apply
glycans_TS_NN <- DA_glycan_apply(contrast.matrix, TIFNIF[,1:Tn], mod_design, 0, 0.05, NULL, FALSE)



# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------






											### VISUALIZATION OF RESULTS ###






# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# GET EXPRESSION MATRICES FOR DA GLYCANS
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

glycans_NC <- c("GP1","GP2","GP3", "GP4","GP10", "GP15", "A1")
glycans_TS <- c("GP1","GP2","GP3", "GP4", "GP7", "GP10", "GP11", "GP15", "A1", "GP0")
glycans_TS_extra <- sort(c(glycans_TS, "G8", "GP9", "GP14", "GP20", "GP23", "coreF"))
glycans_TS_NN <- c("A1", "GP1", "GP2", "GP3", "GP8", "GP9", "GP14","GP18", "GP20", "GP22")

glycans_batch_NC <- batch_corr_NC[rownames(batch_corr_NC) %in% glycans_NC, ]
glycans_raw_NC <- TIFNIF[rownames(TIFNIF) %in% glycans_NC, ]

#glycans_batch_TS <- batch_corr_TS[rownames(batch_corr_TS) %in% glycans_TS, ]
#glycans_raw_TS <- TIFNIF[rownames(TIFNIF) %in% glycans_TS, ]

glycans_batch_TS_extra <- batch_corr_TS[rownames(batch_corr_TS) %in% glycans_TS_extra, ]
glycans_raw_TS_extra <- TIFNIF[rownames(TIFNIF) %in%  glycans_TS_extra, ]

#glycans_TS_NN_raw <- TIFNIF[,1:Tn][rownames(TIFNIF[,1:Tn]) %in% glycans_TS_NN,]  
#glycans_TS_NN_batch <- batch_corr_TS_NN[rownames(batch_corr_TS_NN) %in% glycans_TS_NN,]  



# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# HEATMAP VISUALIZATION
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



# Matrices with colors for Colsidecolors in heatmap
ColorTS <- get_colors(as.character(TS), TS.cols)
ColorTS_NN <- get_colors(as.character(TS_NN), TS.cols[-6])
ColorNC <- get_colors(as.character(NC), NC.cols)
spacer <- as.matrix(rep("white", ncol(TIFNIF)))

# heatmap colors in blue
heat.cols <- colorRampPalette(brewer.pal(9, "YlGnBu"))(n = 300)


# Make pdf with heatmap
pdf("TS.pdf", height = 6, width = 4)
heatmap.plus(as.matrix(scale(glycans_batch_TS_extra , scale = FALSE)), col=heat.cols, Rowv=NULL, hclustfun=function(d) hclust(d, method="ward.D2"), trace="none", labRow=rownames(glycans_batch_TS_extra), labCol=as.character(TS), ColSideColors=cbind(spacer,ColorTS), margins = c(14,8), cexCol=1.2, cexRow = 1.3)
legend(0, 1, legend=c("TS", "normal", "LumA", "LumB", "LumB_HER2_enriched", "HER2", "TNBC", "Apocrine") , fill=c("white", "grey60", "violetred4", "violetred", "orchid2", "purple", "orange", "grey20"), border=FALSE, bty="n", y.intersp = 0.8, cex=1.1)
map <- makecmap(-3:4)
map$colors <- colorRampPalette(brewer.pal(9, "YlGnBu"))(n = length(map$breaks)-1)
hkey(map, x = 0, y = 0, title = "LogFC", stretch = 2)
dev.off()



# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------






							###  TUMOR INFILTRATING LYMPHOCYTE STATUS AND N-GLYCAN ABUNDANCE ###






# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# DATASET - TIF SAMPLES ONLY
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

TIF <- TIFNIF[, 1:Tn]
TIFinfo <- TIFNIFinfo[1:Tn,]

# remove replicate 109.1, as no longer paired analysis
TIF <- TIF[, -which(TIFinfo$ID == "TIF109.1")]
TIFinfo <- TIFinfo[-which(TIFinfo$ID == "TIF109.1"), ]

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# DAA  - TUMOR INFILTRATING LYMPHOCYTES
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# TIL status vector
TILS <- as.factor(as.character(TIFinfo$TILS))

# Simple TIL status, high, low.
TILS <- ifelse(TILS %in% c("T0", "T1"), "low", as.character(TILS))
TILS <- ifelse(TILS %in% c("T2", "T3"), "high", as.character(TILS))

# Remove samples with no TIL info
remove <- c(which(is.na(TIFinfo$TILS)), which(TIFinfo$TILS == "T3_outside_tumor"))

# N-glycans DA between samples with high and low TILs
resTILS <- DA_all_contrasts(TIF, TILS, 0, 0.06, c("high", "low"), remove)



# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# DAA INDIVIDUAL TILS - CD45, CD3, CD4, CD8, CD68
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# remove samples with no CD45 status 
remove <- which(is.na(TIFinfo$CD45))

# CD45 status vector
CD45 <- factor(as.character(TIFinfo$CD45))

resCD45 <- DA_all_contrasts(TIF, CD45, 0, 0.05, c("CD45three", "CD45two", "CD45one"), remove)

# Reform and save data for Figure 2A
#data_ordered <- rbind(resCD45$`GCD45three-GCD45one`[[1]][order(resCD45$`GCD45three-GCD45one`[[1]]$logFC, decreasing=TRUE),], resCD45$`GCD45three-GCD45one`[[2]][order(resCD45$`GCD45three-GCD45one`[[2]]$logFC),])
#Figure2Adata <- data.frame(rownames(data_ordered), data_ordered$logFC, 1/data_ordered$adj.P.Val, log2(1/data_ordered$adj.P.Val))
#colnames(Figure2Adata) <- c("GPs", "logFC", "inverseFDR", "scaledinverseFDR")
#save(Figure2Adata, file="Figure2Adata.Rdata")  




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


# Reform and save data for Figure 2B
#data_ordered <- rbind(resCD4$`GCD4three-GCD4one`[[1]][order(resCD4$`GCD4three-GCD4one`[[1]]$logFC, decreasing=TRUE),], resCD4$`GCD4three-GCD4one`[[2]][order(resCD4$`GCD4three-GCD4one`[[2]]$logFC),])
#Figure2Bdata <- data.frame(rownames(data_ordered), data_ordered$logFC, 1/data_ordered$adj.P.Val, log2(1/data_ordered$adj.P.Val))
#colnames(Figure2Bdata) <- c("GPs", "logFC", "inverseFDR", "scaledinverseFDR")
#save(Figure2Bdata, file="Figure2Bdata.Rdata")  



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

# remove samples with no ER status
remove <- which(is.na(TIFinfo$ER))

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

