# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------




											### LOADING AND PREPARING DATA ###





# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# LOAD DATA AND METADATA
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Glycan data
TIFNIF <-  openxlsx::read.xlsx("TIFNIF_all.xlsx", colNames = TRUE, rowNames = TRUE)

# Metadata
TIFNIFinfo <-  openxlsx::read.xlsx("TIFNIF_all_info.xlsx", colNames = TRUE, rowNames = TRUE)

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
# LOAD SERUM DATA
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

serum <- openxlsx::read.xlsx("serum_glycan_ID.xlsx", colNames = TRUE, rowNames = TRUE)


# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# MATCH TIF WITH SERUMDATA
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


# Get TIF IDs as numbers
s <- as.numeric(unlist(strsplit(TIFinfo$ID, "[^[:digit:]]")))
s <- s[!is.na(s)]
s <- s[-which(s==1)]


# Add serum samples to TIF metadata
TIFinfo$s <- s

# Extract only pairs from metadata with a corresponding serum ID
TIFinfo <- TIFinfo[TIFinfo$s %in% serum$ID, ]

# Extract matching TIF abundance data
TIF <- TIF[,colnames(TIF) %in% TIFinfo$ID]

# Final serum dataset
serum$ID <- NULL
serum <- data.frame(t(serum))


# Log2 tranformation of TIF and serum data
serumlog2 <- log2(serum)
TIFlog2 <- log2(TIF)

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------







											### PEARSONS CORRELATION ANALYSIS ###







# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# PEARSONS CORRELATION ANALYSIS
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Subsets of N-glycans identified as DA in different analysis
glycans_DA <- c("A1", "GP1", "GP2", "GP3", "GP4", "GP10", "GP15")
glycans_DA_all <- c("GP1" , "GP3", "A1",  "GP30",  "GP22",  "GP36",  "GP24",  "GP11",  "G0", "GP40",  "GP39",  "GP17",  "GP43",  "GP10",  "GP23",  "GP15",  "GP9",   "GP4",   "GP8",   "A2",    "GP2",   "GP25", "GP37",  "G2",  "GP26",  "S2",    "GP7",   "coreF", "GP35",  "GP20",  "GP28",  "GP14",  "GP38") 
glycans_DA_TS <- c("GP8", "GP9", "GP14", "GP20", "GP23", "coreF")
glycans_DA_TILS <- c("G0", "S0", "GP5", "S4", "GP41", "GP45")
glycans_DA_CD4CD45 <- c("A4","G1", "G4", "S0", "S4", "GP25", "GP38", "GP37", "GP41", "GP45", "outerarmF")
glycans_DA_survival <- c("GP5", "GP10", "GP23", "GP24", "GP38", "coreF")


# Running function my_correlations on matched TIF and serum using subsets of query N-glycans
DA_top7 <- my_correlation(TIFlog2, serumlog2, glycans_DA, 1)
DA_all <- my_correlation(TIFlog2, serumlog2, glycans_DA_all, 1)
DA_TS <- my_correlation(TIFlog2, serumlog2, glycans_DA_TS, 1)
DA_TILS <- my_correlation(TIFlog2, serumlog2, glycans_DA_TILS, 1)
DA_CD4CD45 <- my_correlation(TIFlog2, serumlog2, glycans_DA_CD4CD45, 1)
DA_survival <- my_correlation(TIFlog2, serumlog2, glycans_DA_survival, 1)

# Reform and save data for Figure 4A and 4B.
#DA_all$GPs <- rownames(DA_all)
#DA_all$Scaled_Inverse_FDR <- log2(1/as.numeric(DA_all$fdr))
#Figure4ABdata <- list(DA_all, TIFlog2, serumlog2)
#save(Figure4ABdata, file="Figure4ABdata.Rdata")




# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------







											### COMPARISON WITH MDG COHORT ###







# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# PREPARE MDG DATASET 
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Load Data
mdg <- read.xlsx("MDG_simple.xlsx", rowNames = TRUE)
# Remove enteries without Subtype
mdg <- mdg[-c(which(is.na(mdg$Subtype))),]

# Vector with Tumor and Normal
mdg.CN <- as.factor(mdg$Tissue)
# Vector with subtypes
mdg.TS <- as.factor(mdg$Subtype)

# Remove metadata to make numeric matrix of abundance values
mdg <- mdg[!colnames(mdg) %in% c("Tissue", "PGR", "ER", "HER2", "Subtype")]
samples <- rownames(mdg)
mdg<- as.matrix(as.data.frame(lapply(mdg, as.numeric)))
mdg <- data.frame(t(mdg))
colnames(mdg) <- samples





# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# DIFFERENTIAL EXPRESSION ANALYSIS WITH MDG TUMOR VS NORMAL
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Rename vector with MDG tumor and normal
NC <- mdg.CN

# Model matrix
mod_design <-  model.matrix(~0+NC)

# Contrast for LIMMA
cotr_N_C <- makeContrasts("NCtumour-NCnormal", levels=mod_design)

# DAA, calling function DA_glycans, logFC cutoff 0 and FDR < 0.05
glycans_MDG_NC <- DA_glycan(cotr_N_C, log2(mdg), mod_design, 0, 0.05)



# Reform and save data for Suplementary Figure 3.
#ColorNC <- get_colors(as.character(NC), c("grey50","red2"))
#spacer <- as.matrix(rep("white", ncol(mdg)))
#Overlap_DA <- mdg[rownames(mdg) %in% c("GP8", "GP9", "GP14", "GP23", "core.F"),]
#Supfig3data <- list(Overlap_DA, cbind(spacer, ColorNC))
#save(Supfig3data, file="Supfig3data.Rdata")




# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# DIFFERENTIAL EXPRESSION ANALYSIS WITH MDG BC SUBTYPES AND NORMAL
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Rename vector with MDG subtypes
TS <- mdg.TS

# Create all combinations of BC subtypes for LIMMA DA analysis
combinations <- data.frame(t(combn(paste0("TS", levels(TS)), 2)))
combinations$contr <- apply(combinations[,colnames(combinations)], 1, paste, collapse = "-")


# Design matrix
mod_design <-  model.matrix(~0+TS)

# Making group contrasts 
contrast.matrix <- eval(as.call(c(as.symbol("makeContrasts"),as.list(as.character(combinations$contr)),levels=list(mod_design))))

# DAA, calling function DA_glycans, logFC cutoff 0 and FDR < 0.05
glycans_MDG_TS <- DA_glycan_apply(contrast.matrix, log2(mdg), mod_design, 0, 0.05, NULL, FALSE)


