# NB! The script N_Glycan_Functions.R must be run in order to make sure any dependencies on packages or custom functions is taking into account. 
# Data for figures are generated from scripts N_Glycan_analysis_paired.R, N_Glycan_Serum_Analysis.R and N_Glycan_Survival_Analysis.R. The code writing out the data is hashtagged out'''

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------





															### FIGURES ###



# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

									### FIGURE 1A ###
									
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

load("Figure1Adata.Rdata")

# Figure 1A
#pdf("Figure1A.pdf", width = 12, height = 6)
#myMDSplot(Figure1Adata[[1]], Figure1Adata[[2]], "", Figure1Adata[[3]])
#dev.off()



# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

									### FIGURE 1B ###

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

load("Figure1Bdata.Rdata")

pdf("Figure1B.pdf", width=10, height = 6)
ggplot(data = Figure1Bdata, aes(x = GPs, y = logFC, fill=scaledinverseFDR)) + geom_bar(stat = 'identity') + scale_fill_gradient2(low="#EDF8B1", mid="#7FCDBB", high="#2C7FB8", space='Lab', name="Inverse FDR") + scale_x_discrete(limits=as.character(Figure1Bdata$GPs)) + theme_bw() + theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_line( size=.1, color="black" )) + theme(axis.text.x = element_text(size=12, color = "black", angle = 90, hjust = 1), axis.text.y = element_text(size=14, color = "black"), axis.title = element_text(size=16, color = "black"), legend.text = element_text(size=16), legend.title = element_text(size=16))
dev.off()





# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

									### FIGURE 2A and FIGURE 2B ###

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

load("Figure2Adata.Rdata")
load("Figure2Bdata.Rdata")


# Barplots TILs
pCD45 <- ggplot(data = Figure2Adata, aes(x = GPs, y = logFC, fill=inverseFDR)) + geom_bar(stat = 'identity') + scale_fill_gradient2(low="white", mid="yellow", high="darkgreen", space='Lab', name="Inverse FDR") + scale_x_discrete(limits=as.character(Figure2Adata$GPs)) + theme_bw() + theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_line( size=.1, color="black" )) + theme(axis.text.x = element_text(size=12, color = "black", angle = 90, hjust = 1), axis.text.y = element_text(size=14, color = "black"), axis.title = element_text(size=14, color = "black"), legend.text = element_text(size=13), legend.title = element_text(size=13))
pCD4 <- ggplot(data = Figure2Bdata, aes(x = GPs, y = logFC, fill=scaledinverseFDR)) + geom_bar(stat = 'identity') + scale_fill_gradient2(low="white", mid="yellow", high="darkgreen", space='Lab', name="Inverse FDR") + scale_x_discrete(limits=as.character(Figure2Bdata$GPs)) + theme_bw() + theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_line( size=.1, color="black" )) + theme(axis.text.x = element_text(size=12, color = "black", angle = 90, hjust = 1), axis.text.y = element_text(size=14, color = "black"), axis.title = element_text(size=14, color = "black"), legend.text = element_text(size=13), legend.title = element_text(size=13))

pdf("Figure2AB.pdf", height = 6, width = 8)
multiplot(pCD45, pCD4, cols = 1)
dev.off()





# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

									### FIGURE 3A ###

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

load("Figure3Adata.Rdata")

pdf("Figure3A.pdf", height = 7, width=12)
my_survival(surv_age_outcome, "")
dev.off()




# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

									### FIGURE 3B ###

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

load("Figure3Bdata_1.Rdata")
load("Figure3Bdata_2.Rdata")


dir.create("survivalcurves")
setwd(paste0(getwd(), "/survivalcurves"))

pdf("Figure3B_1.pdf", height = 7, width = 10)
ggsurvplot(Figure3Bdata_2[[1]], data = Figure3Bdata_1[[1]][1], risk.table = FALSE, font.tickslab = 16, legend.title = "GP5 abundance", legend.labs = c("high", "low"), font.x = 16, font.y = 16, font.legend = 18, xlab = "Years")
dev.off()

pdf("Figure3B_2.pdf", height = 7, width = 10)
ggsurvplot(Figure3Bdata_2[[2]], data = Figure3Bdata_1[[1]][2], risk.table = FALSE, font.tickslab = 16, legend.title = "GP10 abundance", legend.labs = c("high", "low"), font.x = 16, font.y = 16, font.legend = 18, xlab = "Years")
dev.off()

pdf("Figure3B_3.pdf", height = 7, width = 10)
ggsurvplot(Figure3Bdata_2[[3]],  data = Figure3Bdata_1[[1]][3], risk.table = FALSE, font.tickslab = 16, legend.title = "GP23 abundance", legend.labs = c("high", "low"), font.x = 16, font.y = 16, font.legend = 18, xlab = "Years")
dev.off()

pdf("Figure3B_4.pdf", height = 7, width = 10)
ggsurvplot(Figure3Bdata_2[[4]], data = Figure3Bdata_1[[1]][4], risk.table = FALSE, font.tickslab = 16, legend.title = "GP24 abundance", legend.labs = c("high", "low"), font.x = 16, font.y = 16, font.legend = 18, xlab = "Years")
dev.off()

pdf("Figure3B_5.pdf", height = 7, width = 10)
ggsurvplot(Figure3Bdata_2[[5]], data = Figure3Bdata_1[[2]][5], risk.table = FALSE, font.tickslab = 16, legend.title = "GP38 abundance", legend.labs = c("high", "low"), font.x = 16, font.y = 16, font.legend = 18, xlab = "Years")
dev.off()

pdf("Figure3B_6.pdf", height = 7, width = 10)
ggsurvplot(Figure3Bdata_2[[6]], data = Figure3Bdata_1[[2]][6], risk.table = FALSE, font.tickslab = 16, legend.title = "CoreF abundance", legend.labs = c("high", "low"), font.x = 16, font.y = 16, font.legend = 18, xlab = "Years")
dev.off()

setwd("..")




# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

									### FIGURE 4A and 4B ###

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


load("Figure4ABdata.Rdata")

pdf("Figure4A.pdf", height = 8, width = 14)
ggplot(Figure4ABdata[[1]], aes(GPs, cor_coef)) +  geom_point(aes(colour = Scaled_Inverse_FDR), size=7) + scale_colour_gradient(low="lightskyblue1", high="navyblue") + scale_y_continuous(breaks=number_ticks(10)) + theme_bw() +  theme(panel.grid.major.x = element_blank()) + geom_text(aes(label=Figure4ABdata[[1]]$GPs), size=6, hjust = 0.8, vjust=-0.2, color="grey30") + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y = element_text(size=14, color = "black"), axis.title = element_text(size=16, color = "black"), legend.text = element_text(size=16), legend.title = element_text(size=14)) + xlab("") + ylab("Correlation Coefficient") + geom_hline(yintercept=c(0.0, 0.38, -0.45), color=c("grey30","maroon3", "maroon3"))
dev.off()

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

GPs <- c(1,8,9,14,20,23,28,37,38,62)
my_correlation_plots("Figure4B", Figure4ABdata[[2]], Figure4ABdata[[3]], GPs)




# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

									### FIGURE 5 ###

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

load("Figure5data.Rdata")

datf1x <- melt(Figure5data, id.var = 'GPs')
datf1x$value <- as.factor(datf1x$value)

pdf("Figure9.pdf")
ggplot(datf1x, aes(GPs, variable)) + geom_tile(aes(fill = value), colour = "white") + theme_bw()  + scale_fill_manual(values=c("grey80", "cyan3"))+ scale_x_discrete(expand=c(0,0)) + theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) + theme(axis.text.x = element_text(size=13, colour = "black", angle = 90, hjust = 1), axis.text.y = element_text(size=13, colour = "black", angle = 0, hjust = 1), axis.title = element_text(size=13), legend.text = element_text(size=16)) + labs(x = "GP", y = "Trait")
dev.off()





# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

### Supplementary Figure 3 ###

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

load("Supfig3data.Rdata")

# heatmap colors in blue
heat.cols <- colorRampPalette(brewer.pal(9, "YlGnBu"))(n = 300)

pdf("Supfig3.pdf")
heatmap.plus(as.matrix(log2(Supfig3data[[1]])), col=heat.cols, Rowv=NULL, hclustfun=function(d) hclust(d, method="ward.D2"), trace="none", labRow=rownames(Supfig3data[[1]]), labCol="", ColSideColors=Supfig3data[[2]], margins = c(14,8), cexCol=1.2, cexRow = 1.3)
legend(0, 1, legend=c("Normal", "Breast cancer") , fill=c("grey50", "red2"), border=FALSE, bty="n", y.intersp = 0.7, cex=0.9)
map <- makecmap(-3:4)
map$colors <- colorRampPalette(brewer.pal(9, "YlGnBu"))(n = length(map$breaks)-1)
hkey(map, x = 0, y = 0, title = "LogFC", stretch = 2)
dev.off()

