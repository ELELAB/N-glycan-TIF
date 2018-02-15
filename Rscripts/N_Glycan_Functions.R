
# Author: Thilde Bagger Terkelsen
# Contact: thilde@cancer.dk
# Place of employment:Computational Biology Laboratory, Danish Cancer Society Research Center, Copenhagen, Denmark
# Date 29-05-2017

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



													### FUNCTIONS ###






# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# LOAD PACKAGES
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

source("https://bioconductor.org/biocLite.R")
biocLite(c("survcomp", "limma", "sva"))
list.of.packages <- c("limma", "sva", "openxlsx", "ggplot2", "dendextend", "heatmap.plus", "reshape", "gdata", "plyr", "data.table", "RColorBrewer", "squash", "survcomp", "car", "survminer")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, library, character.only=T)


missing.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]


{
if(length(missing.packages) > 0) {print(paste0(missing.packages, " package(s) are missing, install it/these with require(name.of.package) or install.packages(name.of.package). You may need to install one or more packages from bioconductor."))}
stopifnot(length(missing.packages)==0)
}



# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# BOXPLOTS AND VIOLINPLOTS
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Function for generating boxplots: 
# Takes as arguments; 
	# a dataframe
	# a name for the plot (given as a sting)
	# a vector of groups
	# a vector of numbers indicating number of samples in each group
	# a vector of colors for each box (one color for each group)
	
my.boxplots <- function(my.data, my.name, my.group, my.sn, my.cols) {
  png(paste0(my.name,".png"), height = 800, width = 1200)
  p1 <- apply(my.data, 1, function(x) ggplot(,aes(my.group, as.numeric(x))) + geom_boxplot(aes(fill = my.group)) + scale_fill_manual(values=my.cols) + theme_bw() + theme(legend.position="none") + geom_text(data=data.frame(), aes(x=names(c(by(as.numeric(x), my.group, median))), y=c(by(as.numeric(x), my.group, median)), label=my.sn), col='black', size=6) + labs(x = "BC Subtypes", y = "GLycan Abundance"))
  nc <- ceiling(nrow(my.data)/4)
  multiplot(plotlist = p1, cols = nc)
  dev.off()
}

# Function for generating violin plots: 
# Takes as arguments; 
	# a dataframe
	# a name for the plot (given as a sting) 
	# a vector of groups
	# a vector of colors for each box (one color for each group)
	
my.violin <- function(my.data, my.name, my.group, my.cols) {
  pdf(paste0(my.name,".pdf"), height = 6, width = 10)
  p1 <- apply(my.data, 1, function(x) ggplot(,aes(my.group, as.numeric(x))) + geom_violin(aes(fill = my.group), trim = FALSE) + stat_summary(fun.y=median, geom="point", size=2, color="black") + scale_fill_manual(values=my.cols) + theme_bw() + theme(legend.position="none") + theme(axis.text = element_text(size=13, colour = "black"), axis.title = element_text(size=13)) + labs(x = "GP", y = "Log abundance"))
  nc <- ceiling(nrow(my.data)/2)
  multiplot(plotlist = p1, cols = nc)
  dev.off()
}




# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Multidimensional Scaling Plot: 
# Takes as arguments; 
	# a dataframe 
	# a vector of IDs for coloring and labeling (may be the same or different, length should be equal to ncol(dataframe))
	# a vector of colors (one color for each group)
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


myMDSplot <- function(my.data, my.group, my.labels, my.cols) {
  d<-dist(t(my.data))
  fit <- cmdscale(d,eig=TRUE, k=2)
  res<-data.frame(names=rownames(fit$points),M1=fit$points[,1],M2=fit$points[,2])
  p <- ggplot(data=res) 
  p + geom_point(aes(x=M1,y=M2,color=my.group)) + geom_text(aes(x=M1,y=M2, label= my.labels, color=my.group)) + scale_color_manual(values  = my.cols) +
    coord_cartesian(xlim=c(min(res$M1)*1.4,max(res$M1)*1.4)) + theme_bw() + theme(legend.title=element_blank()) + theme(legend.text = element_text(size = 16, face="bold"), axis.title=element_text(size=16,face="bold")) + guides(colour = guide_legend(override.aes = list(size=6))) + theme(axis.text=element_text(size=16, face="bold")) + theme(axis.text = element_text(colour = "black"))
}




# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# FUNCTION TO OBTAIN DIFFERENTIALLY ABUNDANT GLYCANS: 
# Takes as arguments; 
	# a contrast between groups of interest
	# a dataframe, a design matrix with all comparisons
	# cutoffs for logFC and FDR 
	# if blocking than a vector of patient IDs
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


DA_glycan <- function(my.contrast, my.data, my.design, coLFC, coFDR, my.block=NULL) {
  if(is.null(my.block)) {
    fit3 <- eBayes(contrasts.fit(lmFit(my.data, my.design), my.contrast))
  }
  else {
    corfit <- duplicateCorrelation(my.data, my.design, block=my.block) 
    fit3 <- eBayes(contrasts.fit(lmFit(my.data, my.design, block = my.block, correlation=corfit$consensus), my.contrast))
  }
  tt <- toptable(fit3, coef=1, adjust='fdr', number=nrow(my.data))
  
  up <- tt[tt$logFC >= coLFC & tt$adj.P.Val < coFDR, ]
  down <- tt[tt$logFC <= -coLFC & tt$adj.P.Val < coFDR, ]
  
  up$dir <- rep("up", nrow(up))
  down$dir <- rep("down", nrow(down))
  
  final <- list(up, down)
  return(final)
}



# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# FUNCTIONS TO APPLY DIFFERENTIALLY ABUNDANCE ANALYSIS TO ALL COMPARISONS AT ONCE: 
# Takes as arguments; 
	# all contrasts between groups of interest
	# a dataframe
	# a design matrix with all comparisons
	# cutoffs for logFC and FDR
	# if blocking than a vector of patient IDs
	# TRUE/FALSE statment specifying output format, if TRUE the function return a vector of glycan IDs only
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


DA_glycan_apply <- function(my.contrasts, my.data, my.design, coLFC, coFDR, my.block=NULL, my.vector) {
  my.glycans.l <- apply(my.contrasts, 2, function(x) DA_glycan(x, my.data, my.design, coLFC, coFDR, my.block)) 
  if(my.vector == TRUE) {
    my.glycans <- do.call(rbind, lapply(my.glycans.l, function(x) do.call(rbind, x)))
    my.glycans <- unique(do.call(rbind, strsplit(rownames(my.glycans), "[.]"))[,2])
    return(my.glycans)
  }
  else {
    return(my.glycans.l)
  }
}




# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# FUNCTION FOR DA ANALYSIS WITH CLINICAL PARAMETERS. THE FUNCTION CALLS THE "DA_glycan_apply" FROM ABOVE.
# Takes as arguments; 
	# a dataframe
	# a vector of groups do perform contrasts on (same length as ncol(dataframe))
	# a cutoff for logFC and FDR
	# if remove is different from NULL, a vector of indices to remove must be supplied

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

DA_all_contrasts <- function(my.data, my.group, my.logFC, my.FDR, my.levels, my.remove=NULL) {
  if (!is.null(my.remove)) {
    my.data <- my.data[, -my.remove]
    my.group <- my.group[-my.remove]
  }
    G <- factor(as.character(my.group), levels=my.levels)
    combinations<- data.frame(t(combn(paste0("G", levels(G)), 2)))
    combinations$contr <- apply(combinations[,colnames(combinations)], 1, paste, collapse = "-")
    mod_design <-  model.matrix(~0+G)
    contrast.matrix <- eval(as.call(c(as.symbol("makeContrasts"),as.list(as.character(combinations$contr)),levels=list(mod_design))))
    my.DA <- DA_glycan_apply(contrast.matrix, my.data, mod_design, my.logFC, my.FDR, NULL, FALSE)
    return(my.DA)
  }





# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# FUNCTION FOR GETTING HEATMAP COLORS
# Takes as arguments: 
	# a vetcor of groups/labels (a character vector, length of ncol in the matrix to be plotted)
	# a vector with colors to use (a character vector with the length of the number of groups/levels).
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


get_colors <- function(my.truestatus, my.cols) {
  hm_col <- data.frame(levels(as.factor(my.truestatus)), my.cols)
  colnames(hm_col) <- c("status", "mycolor")
  true_status <- data.frame(my.truestatus)
  myorder <- 1:nrow(true_status)
  true_status$order <- myorder
  colnames(true_status) <- c("status", "order")
  col <- merge(true_status, hm_col, by="status", all.x =TRUE)
  col <- col[order(col$order),]
  col$mycolor <- ifelse(is.na(col$mycolor), "black", as.character(col$mycolor)) 
  return(as.matrix(col$mycolor))
}





# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# FUNCTION FOR GENERATING SURVIVAL PLOT
# Takes as arguments;
	# A survival object in the form of a list with a cox regression for each variable (N-glycan).
	# A title of plot	
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

number_ticks <- function(n) {function(limits) pretty(limits, n)}

my_survival <- function(my.survivaldata, my.title) {
  survival_conf <- data.frame(do.call(rbind, lapply(my.survivaldata, function(x) summary(x)$conf.int[nrow(summary(x)$conf.int),c(1,3,4)])))
  colnames(survival_conf) <- c("HR", "lower", "upper")
  survival_pvals <- do.call(rbind, lapply(my.survivaldata, function(x) summary(x)$coefficients[nrow(summary(x)$coefficients),6]))
  survival_fdr <- p.adjust(survival_pvals[,1], method = "fdr", n=nrow(survival_pvals))
  survival_conf$GP <- as.factor(names(my.survivaldata))
  survival_conf$pval <- survival_pvals[,1]
  survival_conf$fdr <- survival_fdr
  survival_conf$sig <- as.factor(ifelse(survival_conf$fdr <=0.05, 1, 0))
  ggplot(survival_conf, aes(x = GP, y=HR)) + geom_point(aes(colour = sig, size = 1/survival_conf$fdr)) + geom_errorbar(aes(ymax = upper, ymin = lower)) + geom_hline(yintercept=1) + scale_x_discrete(limits=as.character(names(my.survivaldata))) + theme_bw() + theme(axis.text.x = element_text(size=13, color = "black", angle = 90, hjust = 1), axis.text.y = element_text(size=12, color = "black"), axis.title = element_text(size=15)) + theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_line( size=.1, color="black" )) + scale_y_continuous(breaks=number_ticks(10)) + ggtitle(my.title) + xlab("GPs") + ylab("Hazard Ratio") + scale_y_continuous(trans = log2_trans(), breaks = trans_breaks("log2", function(x) 2^x),labels = trans_format("log2", math_format(2^.x)))
}





# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# FUNCTION FOR GENERATING SURVIVAL GGPLOT CURVES
# Takes as arguments;
	# a dataframe with abundances
	# a dataframe with survivalinfo
	# indices of features (glycans) to use
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

my_survival_curve <- function(my.data, my.survivaldata, my.index) {
  GPs <- lapply(my.index, function(x) data.frame(as.numeric(my.data[x,]), my.survivaldata$time_to_Outcome_years, my.survivaldata$Outcome_status, my.survivaldata$Age_at_surgery))
  GPs <- lapply(GPs, setNames, c("GP", "years", "status", "age"))
  GPs <- lapply(GPs, function(x) coxph(Surv(years, status) ~ age + GP, data = x))
  return(GPs)
}






# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# FUNCTION FOR CORRELATION ANALYSIS
# Takes as arguments: 
	# two dataframes with values to be correlated. These must have the same dimensions and rownames must the same.
	# list of features to be correlated (e.g. a set of glycans), if all features are to be used then set feature to  rownames(df1) 

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


my_correlation <- function(d1, d2, my.features, my.fdr) {
  
  d1 <- d1[rownames(d1) %in% my.features,]
  d2 <- d2[rownames(d2) %in% my.features,]
  
  my.names <- rownames(d1)
  
  d1 <- as.matrix(as.data.frame(lapply(d1, as.numeric)))
  d2 <- as.matrix(as.data.frame(lapply(d2, as.numeric)))
  
  pear_corr <- sapply(1:nrow(d1), function(i) cor(d1[i,], d2[i,], method = "pearson"))
  pear_corr <- data.frame(pear_corr)
  colnames(pear_corr) <- "cor_coef"
  
  # pearson correlation p-values
  pear_p_val <- sapply(1:nrow(d1), function(i) cor.test(d1[i,], d2[i,], method = "pearson")$p.value)
  pear_p_val <- data.frame(pear_p_val)
  colnames(pear_p_val) <- "pval"
  
  # correction for multiple testing with fdr
  fdr <- data.frame(p.adjust(pear_p_val$pval, method = "fdr"))
  colnames(fdr) <- "fdr"
  
  pear_corr_full <- cbind(pear_corr, pear_p_val, fdr)
  rownames(pear_corr_full) <- my.names
  pear_corr_full <- pear_corr_full[pear_corr_full$fdr <= my.fdr,]
  return(pear_corr_full)
}




# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# FUNCTION FOR CORRELATION SCATTER PLOTS
# Takes as arguments;
	# a name for the plot
	# a dataframe with abundances
	# a dataframe with survivalinfo
	# indices of features (glycans) to use
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

my_correlation_plots <- function(my.name, my.TIFdata, my.serumdata, my.index) {
  GPs <- lapply(my.index, function(x) data.frame(as.numeric(my.TIFdata[x,]), as.numeric(my.serumdata[x,])))
  GPs <- lapply(GPs, setNames, c("TIF", "Serum"))
  p1 <- lapply(GPs, function(x) ggplot(x, aes(TIF, Serum)) + geom_point(shape=1, size=2.5) + theme_bw() + geom_smooth(method=lm) +  theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) + theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12, color = "black"), axis.title = element_text(size=16, color = "black"), legend.text = element_text(size=16)))
  nc <- ceiling(length(GPs)/2)
  pdf(paste0(my.name,".pdf"), height=6, width=12)
  multiplot(plotlist = p1, cols = nc)
  dev.off()
}






# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# FUNCTION FOR MAKING MULTIPLE GGPLOTS IN ONE FIGURE - FUNCTION WAS OBAINED FROM
# http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


