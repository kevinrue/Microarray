# Analysis run to produce the graphs obtained in the draft chapter 1 of my thesis


# Quality control of data for 7 animals (91 samples) ----------------------


# Directory for the analysis of 5 animals
setwd("C:/Users/krue/Documents/Kevin-Logs/20130320_MDM/R_all-91-samples/")

#General purpose library
library(Biobase)

# Load an annotated data frame
targets.91 = read.AnnotatedDataFrame ("phenodata.txt", header = TRUE, row.names = 2, as.is = TRUE)

# Extracts the samples name
SampleNames.91 <- sampleNames(targets.91)

# Directory containing the microarray data files
setwd("C:/Users/krue/Documents/Kevin-Logs/CELfiles/Project/ProjectData/CEL Files detailed")

# Library to manipulate Affymetrix data files
library(affy)

# Imports the raw data
rawdata.91 <- ReadAffy(filenames=pData(targets.91)$FileName, celfile.path=".", phenoData=targets.91, sampleNames=SampleNames.91)

# Directory for the analysis of 5 animals
setwd("C:/Users/krue/Documents/Kevin-Logs/20130320_MDM/R_all-91-samples/")

# Library to import for the quality control step
library(arrayQualityMetrics)

# Quality control
arrayQualityMetrics(expressionset = rawdata.91, outdir = "MDM_raw_log", do.logtransform = TRUE)


# Quality control of data for remaining 5 animals -------------------------


# Directory for the analysis of 5 animals
setwd("C:/Users/krue/Documents/Kevin-Logs/20130320_MDM/R4_5-animals-updated-packages/")

# Load an annotated data frame
targets.65 = read.AnnotatedDataFrame ("phenodata.txt", header = TRUE, row.names = 2, as.is = TRUE)

# Extracts the samples name
SampleNames.65 <- sampleNames(targets.65)

# Directory containing the microarray data files
setwd("C:/Users/krue/Documents/Kevin-Logs/CELfiles/Project/ProjectData/CEL Files detailed")

# Imports the raw data
rawdata.65 <- ReadAffy(filenames=pData(targets.65)$FileName, celfile.path=".", phenoData=targets.65, sampleNames=SampleNames.65)

# Directory for the analysis of 5 animals
setwd("C:/Users/krue/Documents/Kevin-Logs/20130320_MDM/R4_5-animals-updated-packages/") # Directory for the current project

# Quality control of 7 samples
arrayQualityMetrics(expressionset = rawdata, outdir = "MDM_raw_log", do.logtransform = TRUE)


# Normalisation -----------------------------------------------------------


# Library for normalisation
library(farms)

# Normalisation of intensities
farms.65 <- qFarms(rawdata.65)

# One of the two below is required to boxplot FARMS-nomalised expressionSets objects
library(affyPLM)
library(simpleaffy)

# Boxplot of normalised data
boxplot(farms.65, names=SampleNames.65, col = c(rep("pink",65)),las=2, cex.axis=0.53)

# Histogram and legend of normalised data
hist(farms.65, col = length(SampleNames.65):1, lty=length(SampleNames.65):1)
legend("topright", legend=row.names(pData(farms.65)), col=length(SampleNames.65):1, lty=length(SampleNames.65):1, cex=0.60)


# Informative probe sets --------------------------------------------------


# Find informative genes from the farms normalised probe sets using INIcalls
INIs.65 <- INIcalls(farms.65) 
# List of informative probesets, as another exprSet
farms_informative.65 <- getI_Eset(INIs.65) 
# Cleanup
rm(INIs.65)


# Differential gene expression --------------------------------------------


# Creates the design matrix (package stats)
design.65.DGE <- model.matrix(~0 + factor(targets.65$Group, levels = unique(targets.65$Group)))

# Renames the columns of the design matrix
colnames(design.65.DGE) <- unique(targets.65$Group)

# Design the matrix listing the contrasts to compute
contrasts.DGE.names <- c("BOVIS_2HR-MAP_2HR","BOVIS_6HR-MAP_6HR","BOVIS_24HR-MAP_24HR",
            "BOVIS_2HR-BCG_2HR","BOVIS_6HR-BCG_6HR","BOVIS_24HR-BCG_24HR",
            "BOVIS_2HR-CN_2HR","BOVIS_6HR-CN_6HR","BOVIS_24HR-CN_24HR",
            "MAP_2HR-BCG_2HR","MAP_6HR-BCG_6HR","MAP_24HR-BCG_24HR",
            "MAP_2HR-CN_2HR","MAP_6HR-CN_6HR","MAP_24HR-CN_24HR",
            "BCG_2HR-CN_2HR","BCG_6HR-CN_6HR","BCG_24HR-CN_24HR")

# Loads library
library(limma)

# Make contrast matrix (limma)
contrasts.DGE <- makeContrasts(contrasts=contrasts.DGE.names, levels=design.65.DGE) 

# Required step for paired analysis (limma, statmod)
# Note: based on entire FARMS data > informative probe sets
array_correlation <- duplicateCorrelation(farms.65, design=design.65.DGE, ndups=1, block=targets.65$Animal)

# Fit a linear model to the FARMS normalized data, call it farms_fit 
farms_fit = lmFit(farms.65, design = design.65.DGE, ndups=1, cor=array_correlation$consensus, block=targets.65$Animal)

#Apply contrasts to the linear model farms_fit
farms_fit_contrast <- contrasts.fit(farms_fit, contrasts.DGE)

# Filter farms_fit_contrast using the informative genes found earlier
farms.contrast.informative <- farms_fit_contrast[rownames(farms_fit_contrast) %in% rownames(exprs(farms_informative.65)),] 

# Below: Given a series of related parameter estimates and standard errors, compute moderated t-statistics, moderated F-statistic,
# and log-odds of differential expression by empirical Bayes shrinkage of the standard errors towards a common value.
ebayes.informative <- eBayes(farms.contrast.informative)

# Write DEgenes to an outfile
setwd(dir="C:/Users/krue/Documents/Kevin-Logs/20130320_MDM/Writing thesis chapter/")
save(ebayes.informative, file = "ebayes.informative.RData") 

# differential expression information of informative probes stored ina file
nrow(ebayes.informative)

# Cleanup
rm(design.65.DGE, contrasts.DGE, array_correlation, farms_fit, farms_fit_contrast, farms.contrast.informative)


# Number of genes up and down regulated (Including function) --------------


# Function returning count of ups and downs in dataset
UpsAndDowns <- function(ebayes, contrasts)
{
  library(package=limma, quietly=T)
  up.down.count = data.frame()
  for(contrast in contrasts){
    DE = topTable(ebayes, coef= contrast, number=nrow(ebayes.informative), adjust.method="BH")
    # Count genes passing thresholds on log-fold change and p-value
    up.FC0p05 = sum(DE$logFC > 0 & DE$adj.P.Val < 0.05)
    down.FC0p05 = sum(DE$logFC < 0 & DE$adj.P.Val < 0.05)
    DE.FC0p05 = sum(DE$adj.P.Val < 0.05)
    nDE.FC0p05 = sum(DE$adj.P.Val > 0.05)
    # put data in data frame 
    up.down.count = rbind(up.down.count,
                          list(Up=up.FC0p05,Down=down.FC0p05,
                               Total.DE=DE.FC0p05,Total.nDE=nDE.FC0p05))
  }
  # renames rows with sample name
  rownames(up.down.count) = contrasts
  # returns the dataset
  return(up.down.count)
}
# Counts of ups and downs in each contrast
ups.and.downs.counts = UpsAndDowns(ebayes=ebayes.informative, contrasts=contrasts.DGE.names)
# Ups and down for each contrast directly written in file
write.table(x=ups.and.downs.counts, file="UpsAndDowns.txt", quote=F, sep="\t", row.names=T, col.names=T)
# Cleanup
rm(ups.and.downs.counts)


# Boolean data frame of whether probes DE in contrasts --------


# Function returning whether each probe is DE in each contrast
DEprobesets <- function(ebayes, contrasts, probesets)
{
  DE.table = data.frame(ID=probesets)
  for(contrast in contrasts){
    topT = topTable(ebayes, coef= contrast, number=nrow(ebayes.informative), adjust.method="BH")
    # Checks whether each genes is DE for that contrast
    topT$bool = topT$adj.P.Val < 0.05
    # Renames the column with the contrast name
    names(topT)[names(topT)=="bool"] = contrast
    # put data in data frame
    DE.table = merge(x=DE.table, y=topT[c("ID",contrast)], by="ID", sort=TRUE)
  }
  return(DE.table)
}
# Return whether each probe is DE in each contrast
DE.probe.contrast = DEprobesets(ebayes=ebayes.informative,contrasts=contrasts.DGE.names,
                                probesets=rownames(assayData(farms_informative.65)$exprs))
# Cleanup
rm(DE.probe.contrast)


# Venn diagrams (limma: boring B&W)  ---------------------------------

# Venn diagram of all mycobacteria-CN at 2 hpi
vennDiagram(object=DE.probe.contrast[,c("MAP_2HR-CN_2HR","BCG_2HR-CN_2HR","BOVIS_2HR-CN_2HR")],
            fill = c("red","green","blue"))
# Venn diagram of all mycobacteria-CN at 6 hpi
vennDiagram(object=DE.probe.contrast[,c("MAP_6HR-CN_6HR","BCG_6HR-CN_6HR","BOVIS_6HR-CN_6HR")])
# Venn diagram of all mycobacteria-CN at 24 hpi
vennDiagram(object=DE.probe.contrast[,c("MAP_24HR-CN_24HR","BCG_24HR-CN_24HR","BOVIS_24HR-CN_24HR")])


# List of significant probe sets in each contrast  ------------

library(limma)
#
decide.tests = decideTests(object=ebayes.informative)
#
summary(decide.tests)
head(decide.tests)
#
significant.probesets = list(
  "MAP_2HR-CN_2HR"=rownames(decide.tests[decide.tests[,"MAP_2HR-CN_2HR"] != 0,]),
  "BCG_2HR-CN_2HR"=rownames(decide.tests[decide.tests[,"BCG_2HR-CN_2HR"] != 0,]),
  "BOVIS_2HR-CN_2HR"=rownames(decide.tests[decide.tests[,"BOVIS_2HR-CN_2HR"] != 0,]),
  "MAP_2HR-BCG_2HR"=rownames(decide.tests[decide.tests[,"MAP_2HR-BCG_2HR"] != 0,]),
  "BOVIS_2HR-BCG_2HR"=rownames(decide.tests[decide.tests[,"BOVIS_2HR-BCG_2HR"] != 0,]),
  "BOVIS_2HR-MAP_2HR"=rownames(decide.tests[decide.tests[,"BOVIS_2HR-MAP_2HR"] != 0,]),
  "MAP_6HR-CN_6HR"=rownames(decide.tests[decide.tests[,"MAP_6HR-CN_6HR"] != 0,]),
  "BCG_6HR-CN_6HR"=rownames(decide.tests[decide.tests[,"BCG_6HR-CN_6HR"] != 0,]),
  "BOVIS_6HR-CN_6HR"=rownames(decide.tests[decide.tests[,"BOVIS_6HR-CN_6HR"] != 0,]),
  "MAP_6HR-BCG_6HR"=rownames(decide.tests[decide.tests[,"MAP_6HR-BCG_6HR"] != 0,]),
  "BOVIS_6HR-BCG_6HR"=rownames(decide.tests[decide.tests[,"BOVIS_6HR-BCG_6HR"] != 0,]),
  "BOVIS_6HR-MAP_6HR"=rownames(decide.tests[decide.tests[,"BOVIS_6HR-MAP_6HR"] != 0,]),
  "MAP_24HR-CN_24HR"=rownames(decide.tests[decide.tests[,"MAP_24HR-CN_24HR"] != 0,]),
  "BCG_24HR-CN_24HR"=rownames(decide.tests[decide.tests[,"BCG_24HR-CN_24HR"] != 0,]),
  "BOVIS_24HR-CN_24HR"=rownames(decide.tests[decide.tests[,"BOVIS_24HR-CN_24HR"] != 0,]),
  "MAP_24HR-BCG_24HR"=rownames(decide.tests[decide.tests[,"MAP_24HR-BCG_24HR"] != 0,]),
  "BOVIS_24HR-BCG_24HR"=rownames(decide.tests[decide.tests[,"BOVIS_24HR-BCG_24HR"] != 0,]),
  "BOVIS_24HR-MAP_24HR"=rownames(decide.tests[decide.tests[,"BOVIS_24HR-MAP_24HR"] != 0,])
)
# Cleanup
rm(decide.tests)
rm(significant.probesets)


# Venn diagrams (Vennerable: fails for 24 h) ----------------------------------------------
# Note that overlap are so small the layout is clumsy
# 6HR timepoint
plot(Venn(Sets=significant.probesets[c("MAP_2HR-CN_2HR",
                                       "BCG_2HR-CN_2HR",                                       
                                       "BOVIS_2HR-CN_2HR")]))
# 6HR timepoint
plot(Venn(Sets=significant.probesets[c("BCG_6HR-CN_6HR",
                                       "MAP_6HR-CN_6HR",
                                       "BOVIS_6HR-CN_6HR")]))

# 24HR timepoint
## Fails (likely due to too many zeros)
# plot(Venn(Sets=significant.probesets[c("BCG_24HR-CN_24HR",
#                                        "MAP_24HR-CN_24HR",
#                                        "BOVIS_24HR-CN_24HR")]))


# Venn diagrams (VennDiagram: good!) ---------------------------------------------
# Loads library
library(VennDiagram)
# 2 hours (written in file)
venn.diagram(cex=2, cat.just=list(c(0.35,-0.4),c(0.65,-0.4),c(0.5,0)),
             fill = c("green", "yellow", "darkorchid1"),
             x=significant.probesets[c("MAP_2HR-CN_2HR",
                                     "BCG_2HR-CN_2HR",
                                     "BOVIS_2HR-CN_2HR")],
             filename="VennDiagram.mycobacteria_control_2h.tiff")
# 6 hours (written in file)
venn.diagram(cex=2, cat.just=list(c(0.35,-0.4),c(0.65,-0.4),c(0.5,0)),
             fill = c("green", "yellow", "darkorchid1"),
             x=significant.probesets[c("MAP_6HR-CN_6HR",
                                       "BCG_6HR-CN_6HR",
                                       "BOVIS_6HR-CN_6HR")],
             filename="VennDiagram.mycobacteria_control_6h.tiff")
# 24 hours (written in file)
venn.diagram(cex=c(2,1,1,1,0.8,1,1), cat.pos=c(0,0,270), cat.dist=c(0.02,0.015,0.063),
             fill = c("darkorchid1", "yellow", "green"),
             x=significant.probesets[c("BOVIS_24HR-CN_24HR",
                                       "BCG_24HR-CN_24HR",
                                       "MAP_24HR-CN_24HR")],
             filename="VennDiagram.mycobacteria_control_24h.tiff")


# Data frame of fold changes for each probe in each contrast --------------


library(limma) # topTable()
# Function returning the fold change in each contrast
extractFoldChanges <- function(ebayes, contrasts, probesets)
{
  FC.table = data.frame(ID=probesets)
  for(contrast in contrasts){
    topT = topTable(ebayes, coef= contrast, number=nrow(ebayes.informative), adjust.method="BH")
    # Renames the column with the contrast name
    names(topT)[names(topT)=="logFC"] = contrast
    # put data in data frame
    FC.table = merge(x=FC.table, y=topT[c("ID",contrast)], by="ID", sort=TRUE)
  }
  return(FC.table)
}
# Obtains the fold change values for each contrast
FC.probe.contrast = extractFoldChanges(ebayes=ebayes.informative,contrasts=contrasts.DGE.names,
                                       probesets=rownames(assayData(farms_informative.65)$exprs))
# Cleanup
rm(FC.probe.contrast)


# Multiple correlation of log2 fold changes ----------------------------------------


# Appropriate library
library(PerformanceAnalytics)
# Multiple pairwise correlation plots and R-values
chart.Correlation(R=FC.probe.contrast[,
      c("MAP_2HR-CN_2HR","MAP_6HR-CN_6HR", "MAP_24HR-CN_24HR",
        "BCG_2HR-CN_2HR","BCG_6HR-CN_6HR","BCG_24HR-CN_24HR",
        "BOVIS_2HR-CN_2HR","BOVIS_6HR-CN_6HR","BOVIS_24HR-CN_24HR")])


# Separate correlation plots of log2 fold changes ----------------------------------------


# Zoomed correlation graphs
library(ggplot2)
# 2 hours BOV/MAP
ggplot.graph <- qplot(FC.probe.contrast$"BOVIS_2HR-CN_2HR",
                      FC.probe.contrast$"MAP_2HR-CN_2HR",
           xlab = "BOVIS_2HR-CN_2HR", ylab="MAP_2HR-CN_2HR",
           xlim=c(-3, 6), ylim=c(-3, 6))
ggplot.graph
ggplot.graph + stat_smooth(method="lm", se=TRUE)
# 2 hours BOV/BCG
ggplot.graph <- qplot(FC.probe.contrast$"BOVIS_2HR-CN_2HR",
                      FC.probe.contrast$"BCG_2HR-CN_2HR",
                      xlab = "BOVIS_2HR-CN_2HR", ylab="BCG_2HR-CN_2HR",
                      xlim=c(-3, 6), ylim=c(-3, 6))
# Initialises the output file
Cairo(1000, 1000, file="FC.correlation/BOVIS.BCG.FC.correlation.6hpi.png", type="png", bg="white")
# Produces the plot
ggplot.graph + stat_smooth(method="lm", se=TRUE) + theme(axis.title.x = element_text(size=20), axis.text.x  = element_text(face="bold", size=16),
                                                         axis.title.y = element_text(size=20), axis.text.y  = element_text(face="bold", size=16))
# Saves teh output file
dev.off()


# 2 hours MAP/BCG
ggplot.graph <- qplot(FC.probe.contrast$"MAP_2HR-CN_2HR",
                      FC.probe.contrast$"BCG_2HR-CN_2HR",
                      xlab = "MAP_2HR-CN_2HR", ylab="BCG_2HR-CN_2HR",
                      xlim=c(-3, 6), ylim=c(-3, 6))
ggplot.graph + stat_smooth(method="lm", se=TRUE)
# 6 hours BOV/MAP
ggplot.graph <- qplot(FC.probe.contrast$"BOVIS_6HR-CN_6HR",
                      FC.probe.contrast$"MAP_6HR-CN_6HR",
                      xlab = "BOVIS_6HR-CN_6HR", ylab="MAP_6HR-CN_6HR",
                      xlim=c(-3, 6), ylim=c(-3, 6))
ggplot.graph + stat_smooth(method="lm", se=TRUE)
# 6 hours BOV/BCG
ggplot.graph <- qplot(FC.probe.contrast$"BOVIS_6HR-CN_6HR",
                      FC.probe.contrast$"BCG_6HR-CN_6HR",
                      xlab = "BOVIS_6HR-CN_6HR", ylab="BCG_6HR-CN_6HR",
                      xlim=c(-3, 6), ylim=c(-3, 6))
ggplot.graph + stat_smooth(method="lm", se=TRUE) + theme(axis.title.x = element_text(size=20), axis.text.x  = element_text(face="bold", size=16),
                                                         axis.title.y = element_text(size=20), axis.text.y  = element_text(face="bold", size=16))
# 6 hours MAP/BCG
ggplot.graph <- qplot(FC.probe.contrast$"MAP_6HR-CN_6HR",
                      FC.probe.contrast$"BCG_6HR-CN_6HR",
                      xlab = "MAP_6HR-CN_6HR", ylab="BCG_6HR-CN_6HR",
                      xlim=c(-3, 6), ylim=c(-3, 6))
ggplot.graph + stat_smooth(method="lm", se=TRUE)+ theme(axis.title.x = element_text(size=20), axis.text.x  = element_text(face="bold", size=16),
                                                        axis.title.y = element_text(size=20), axis.text.y  = element_text(face="bold", size=16))
# 24 hours BOV/MAP
ggplot.graph <- qplot(FC.probe.contrast$"BOVIS_24HR-CN_24HR",
                      FC.probe.contrast$"MAP_24HR-CN_24HR",
                      xlab = "BOVIS_24HR-CN_24HR", ylab="MAP_24HR-CN_24HR",
                      xlim=c(-3, 6), ylim=c(-3, 6))
ggplot.graph + stat_smooth(method="lm", se=TRUE)
# 24 hours BOV/BCG
ggplot.graph <- qplot(FC.probe.contrast$"BOVIS_24HR-CB_24HR",
                      FC.probe.contrast$"BCG_24HR-CN_24HR",
                      xlab = "BOVIS_24HR-CN_24HR", ylab="BCG_24HR-CN_24HR",
                      xlim=c(-3, 6), ylim=c(-3, 6))
ggplot.graph + stat_smooth(method="lm", se=TRUE)
# 24 hours MAP/BCG
ggplot.graph <- qplot(FC.probe.contrast$"MAP_24HR-CN_24HR",
                      FC.probe.contrast$"BCG_24HR-CN_24HR",
                      xlab = "MAP_24HR-CN_24HR", ylab="BCG_24HR-CN_24HR",
                      xlim=c(-3, 6), ylim=c(-3, 6))
ggplot.graph + stat_smooth(method="lm", se=TRUE)


# Get intercept and slope of linear fit
FC.correlations.lm = data.frame()
FC.correlations.lm = rbind(FC.correlations.lm, coef(lm(`BOVIS_2HR-CN_2HR`~`MAP_2HR-CN_2HR`, data = FC.probe.contrast)))
FC.correlations.lm = rbind(FC.correlations.lm, coef(lm(`BOVIS_6HR-CN_6HR`~`MAP_6HR-CN_6HR`, data = FC.probe.contrast)))
FC.correlations.lm = rbind(FC.correlations.lm, coef(lm(`BOVIS_24HR-CN_24HR`~`MAP_24HR-CN_24HR`, data = FC.probe.contrast)))
FC.correlations.lm = rbind(FC.correlations.lm, coef(lm(`BOVIS_2HR-CN_2HR`~`BCG_2HR-CN_2HR`, data = FC.probe.contrast)))
FC.correlations.lm = rbind(FC.correlations.lm, coef(lm(`BOVIS_6HR-CN_6HR`~`BCG_6HR-CN_6HR`, data = FC.probe.contrast)))
FC.correlations.lm = rbind(FC.correlations.lm, coef(lm(`BOVIS_24HR-CN_24HR`~`BCG_24HR-CN_24HR`, data = FC.probe.contrast)))
FC.correlations.lm = rbind(FC.correlations.lm, coef(lm(`MAP_2HR-CN_2HR`~`BCG_2HR-CN_2HR`, data = FC.probe.contrast)))
FC.correlations.lm = rbind(FC.correlations.lm, coef(lm(`MAP_6HR-CN_6HR`~`BCG_6HR-CN_6HR`, data = FC.probe.contrast)))
FC.correlations.lm = rbind(FC.correlations.lm, coef(lm(`MAP_24HR-CN_24HR`~`BCG_6HR-CN_6HR`, data = FC.probe.contrast)))
# Renames the rows in the same order
row.names(FC.correlations.lm) =c("BOVIS_2HR-CN_2HR~MAP_2HR-CN_2HR","BOVIS_6HR-CN_6HR~MAP_6HR-CN_6HR","BOVIS_24HR-CN_24HR~MAP_24HR-CN_24HR",
                                 "BOVIS_2HR-CN_2HR~BCG_2HR-CN_2HR","BOVIS_6HR-CN_6HR~BCG_6HR-CN_6HR","BOVIS_24HR-CN_24HR~BCG_24HR-CN_24HR",
                                 "MAP_2HR-CN_2HR~BCG_2HR-CN_2HR","MAP_6HR-CN_6HR~BCG_6HR-CN_6HR","MAP_24HR-CN_24HR~BCG_24HR-CN_24HR")

# Renames the columns
colnames(FC.correlations.lm) = c("Intercept","Slope")
# Writes the table in a file
write.table(x=FC.correlations.lm, file="FC.correlations.lm.txt", quote=F, sep="\t", dec=".", row.names=T, col.names=T)

# Cleanup
rm(ggplot.graph, FC.correlations.lm)


# Principal components analysis -------------------------------------------


# PCA calculation (fast)
pca.65 <- prcomp(data.frame(t(assayData(farms_informative.65)$exprs)),  scale=T)
# Cleanup line
rm(pca.65)


# PCA coloring code -------------------------------------------------------


# 13 colors to use for the different TimePoint-Treatment groups
time.treatment.col = c("purple","grey80", "grey40", "black",
                          "#FFAAAA", "#FF5555", "#FF0000",
                          "olivedrab1", "olivedrab3", "olivedrab4",
                          "#BDC9EE","#7B94DD","#3A5FCD") # colors for the different groups
time.pch = c(1,17,15,18)
# 13 colors to use for the different TimePoint-Treatment groups
time.cex = c(1.5,1.5,1.5,2)
# Cleanup line
rm(time.treatment.col, time.pch, time.cex)


# PC1-PC2 of FARMS normalised intensities --------


# Saves the default graphical parameters for later restoration
opar = par()

# Plot of PC1 and PC2 using FARMS normalised informative probe sets
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE) # Add extra space to right of plot area; change clipping to figure
plot(pca.65$x,
     cex=time.cex[as.numeric(
       factor(pData(farms_informative.65)$TimePoint,
              levels=unique(pData(farms_informative.65)$TimePoint),
              ordered=T))],
     # Shapes: 0HR-circle, 2HR-triangle, 6HR-square, 24HR-diamond
     pch=time.pch[as.numeric(
       factor(pData(farms_informative.65)$TimePoint,
              levels=unique(pData(farms_informative.65)$TimePoint),
              ordered=T))],
     # Color of data points
     col = time.treatment.col[as.numeric(
               factor(x=pData(farms_informative.65)$Group,
                      levels=unique(pData(farms_informative.65)$Group),
                      ordered=T))]
)
library(calibrate)
textxy(pca.65$x[,1], pca.65$x[,2], pData(farms_informative.65)$Animal,cx=.8)
legend("right",
       inset = c(-0.14,0),
       legend = unique(pData(farms_informative.65)$Group),
       pch = time.pch[as.numeric(factor(x=aggregate(formula=TimePoint~Group,data=pData(farms_informative.65), FUN=function(x) x[1])$TimePoint[unique(
         tapply(X=pData(farms_informative.65)$TimePoint,INDEX=as.factor(pData(farms_informative.65)$Group)))],
                             levels=unique(pData(farms_informative.65)$TimePoint),
                             ordered=T))],
       col = time.treatment.col,
       cex = 1) # Add legend to top right, outside plot region

# Restore original graphical parameters
par = opar
# Cleanup
rm(par, opar)


# PC3-PC4 of FARMS normalised intensities --------


# Saves the default graphical parameters for later restoration
opar = par()

# Plot of PC1 and PC2 using FARMS normalised informative probe sets
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE) # Add extra space to right of plot area; change clipping to figure
plot(x=pca.65$x[,3], y=pca.65$x[,4],
     xlab="PC3",ylab="PC4",
     cex=time.cex[as.numeric(
       factor(pData(farms_informative.65)$TimePoint,
              levels=unique(pData(farms_informative.65)$TimePoint),
              ordered=T))],
     # Shapes: 0HR-circle, 2HR-triangle, 6HR-square, 24HR-diamond
     pch=time.pch[as.numeric(
       factor(pData(farms_informative.65)$TimePoint,
              levels=unique(pData(farms_informative.65)$TimePoint),
              ordered=T))],
     # Color of data points
     col = time.treatment.col[as.numeric(
       factor(x=pData(farms_informative.65)$Group,
              levels=unique(pData(farms_informative.65)$Group),
              ordered=T))]
)
library(calibrate)
#
textxy(pca.65$x[,3], pca.65$x[,4], pData(farms_informative.65)$Animal,cx=.8)
legend("right",
       inset = c(-0.14,0),
       legend = unique(pData(farms_informative.65)$Group),
       pch = time.pch[as.numeric(factor(x=aggregate(formula=TimePoint~Group,data=pData(farms_informative.65), FUN=function(x) x[1])$TimePoint[unique(
         tapply(X=pData(farms_informative.65)$TimePoint,INDEX=as.factor(pData(farms_informative.65)$Group)))],
                                        levels=unique(pData(farms_informative.65)$TimePoint),
                                        ordered=T))],
       col = time.treatment.col,
       cex = 0.9) # Add legend to top right, outside plot region

# Restore original graphical parameters
par = opar
# Cleanup
rm(par, opar)

# Percentage of variance explained by each principal component
summary(pca.65)


# PC5-PC6 of FARMS normalised intensities --------


# Saves the default graphical parameters for later restoration
opar = par()

# Plot of PC1 and PC2 using FARMS normalised informative probe sets
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE) # Add extra space to right of plot area; change clipping to figure
plot(x=pca.65$x[,5], y=pca.65$x[,6],
     xlab="PC5",ylab="PC6",
     cex=time.cex[as.numeric(
       factor(pData(farms_informative.65)$TimePoint,
              levels=unique(pData(farms_informative.65)$TimePoint),
              ordered=T))],
     # Shapes: 0HR-circle, 2HR-triangle, 6HR-square, 24HR-diamond
     pch=time.pch[as.numeric(
       factor(pData(farms_informative.65)$TimePoint,
              levels=unique(pData(farms_informative.65)$TimePoint),
              ordered=T))],
     # Color of data points
     col = time.treatment.col[as.numeric(
       factor(x=pData(farms_informative.65)$Group,
              levels=unique(pData(farms_informative.65)$Group),
              ordered=T))]
)
library(calibrate)
#
textxy(pca.65$x[,5], pca.65$x[,6], pData(farms_informative.65)$Animal,cx=.8)
legend("right",
       inset = c(-0.14,0),
       legend = unique(pData(farms_informative.65)$Group),
       pch = time.pch[as.numeric(factor(x=aggregate(formula=TimePoint~Group,data=pData(farms_informative.65), FUN=function(x) x[1])$TimePoint[unique(
         tapply(X=pData(farms_informative.65)$TimePoint,INDEX=as.factor(pData(farms_informative.65)$Group)))],
                                        levels=unique(pData(farms_informative.65)$TimePoint),
                                        ordered=T))],
       col = time.treatment.col,
       cex = 0.9) # Add legend to top right, outside plot region

# Restore original graphical parameters
par = opar
# Cleanup
rm(par, opar)

# Percentage of variance explained by each principal component
summary(pca.65)


# PCA no time zero --------------------------------------------------------

farms_informative.65_no0H = farms_informative.65[,-grep(pattern="0H", x=rownames(pData(targets.65)))]
targets.65_no0H = targets.65[-grep(pattern="0H", x=rownames(pData(targets.65))),]
pca.65.no0H = prcomp(data.frame(t(assayData(farms_informative.65_no0H)$exprs)),  scale=T)


# PC1-PC2 of FARMS normalised intensities no time 0H --------


# Saves the default graphical parameters for later restoration
opar = par()

# Plot of PC1 and PC2 using FARMS normalised informative probe sets
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE) # Add extra space to right of plot area; change clipping to figure
plot(pca.65.no0H$x,
     cex=time.cex[as.numeric(
       factor(pData(farms_informative.65_no0H)$TimePoint,
              levels=unique(pData(farms_informative.65)$TimePoint),
              ordered=T))],
     # Shapes: 0HR-circle, 2HR-triangle, 6HR-square, 24HR-diamond
     pch=time.pch[as.numeric(
       factor(pData(farms_informative.65_no0H)$TimePoint,
              levels=unique(pData(farms_informative.65)$TimePoint),
              ordered=T))],
     # Color of data points
     col = time.treatment.col[as.numeric(
       factor(x=pData(farms_informative.65_no0H)$Group,
              levels=unique(pData(farms_informative.65)$Group),
              ordered=T))]
)
library(calibrate)
textxy(pca.65.no0H$x[,1], pca.65.no0H$x[,2], pData(farms_informative.65)$Animal,cx=.8)
legend("right",
       inset = c(-0.14,0),
       legend = unique(pData(farms_informative.65_no0H)$Group),
       pch = time.pch[1+as.numeric(factor(x=aggregate(formula=TimePoint~Group,data=pData(farms_informative.65_no0H), FUN=function(x) x[1])$TimePoint[unique(
         tapply(X=pData(farms_informative.65_no0H)$TimePoint,INDEX=as.factor(pData(farms_informative.65_no0H)$Group)))],
         levels=unique(pData(farms_informative.65_no0H)$TimePoint),
         ordered=T))],
       col = time.treatment.col[-c(1)],
       cex = 1) # Add legend to top right, outside plot region

# Restore original graphical parameters
par = opar
# Cleanup
rm(par, opar)


# Annotation of probe sets with cross-reference IDs -----------------------


AnnotateProbesTable = function(df, annotPkg="bovine.db", probeCol="ID",
                               gene.symbol=TRUE, ENSEMBL=TRUE, ENTREZ=TRUE)
{
  # Annotates a data frame containing probe IDs with Official gene symbol,
  # ENTREZ gene ID, and ENSEMBL gene ID.
  #
  # Args:
  #   df: Data frame containing the probe ID.
  #   annotPkg: Annotation package to use. Default is "bovine.db"
  #   probeCol: identifier of the column containing the probe IDs. Default is "ID".
  #   symbol: If TRUE, annotates with official gene symbol.
  #   ENSEMBL: If TRUE, annotates with ENSEMBL gene ID.
  #   ENTREZ: If TRUE, annotates with ENTREZ gene ID.
  #
  # Returns:
  #   The input data.frame with annotation columns appended.
  #
  n = nrow(df)
  df$rankFC <- 1:n # Keeps track of original row order
  # Packages required to annotate
  library(annotate, quietly=TRUE)
  # Fancy way to pass the library to load as an argument
  library(annotPkg, character.only=TRUE, quietly=TRUE)
  # Error handling
  if (!probeCol %in% colnames(df)){
    stop("Column does not exist in data.frame: ",
         probeCol, ".")
  }
  # Obtains the gene symbol annotations
  if (gene.symbol){
    df = cbind(df, gene.symbol=getSYMBOL(df[,probeCol], annotPkg))
  }
  # Obtains the ENSEMBL annotations
  if (ENSEMBL){
    df = merge(df, toTable(bovineENSEMBL[df[,probeCol]]), by.x=probeCol, by.y="probe_id",
               all.x=TRUE)
  }
  # Obtains the ENTREZ annotations
  if (ENTREZ){
    df = merge(df, toTable(bovineENTREZID[df[,probeCol]]), by.x=probeCol, by.y="probe_id",
               all.x=TRUE)
  }
  # Restores the original row order and rownames (messed up by "merge")
  df = df[order(df$rankFC),]
  # Discards the ranking column
  df$rankFC = NULL
  return(df)
}


# F-value of animal effect ------------------------------------------------


# Obtain the F-value for each probe set (2 minutes)
one.way.anova.animal.effect = data.frame(F=apply(X=assayData(farms_informative.65)$exprs, MARGIN=1, FUN=function(x){
                         oneway.test(formula=expr~animal, data=cbind(expr=x, animal=as.factor(targets.65$Animal)))$statistic}))
# Obtain the p-value for each probe set (2 minutes)
one.way.anova.animal.effect$p.value = apply(X=assayData(farms_informative.65)$exprs, MARGIN=1, FUN=function(x){
                                              oneway.test(formula=expr~animal, data=cbind(expr=x, animal=as.factor(targets.65$Animal)))$p.value})
one.way.anova.animal.effect = one.way.anova.animal.effect[order(one.way.anova.animal.effect$F, decreasing=T),]
# For annotation purpose, a column must contain the probe set ID
one.way.anova.animal.effect$ID = rownames(one.way.anova.animal.effect)
# Annotate with Gene Symbol (ENTREZ and ENSEMBL IDs can be obtained too)
one.way.anova.animal.effect.GeneSymbol = AnnotateProbesTable(df=one.way.anova.animal.effect, ENSEMBL=F, ENTREZ=F)
# Moves to folder dedicated to this analysis
setwd(dir="C:/Users/krue/Documents/Kevin-Logs/20130320_MDM/Writing thesis chapter/")
# Write this document in a 
write.table(x=one.way.anova.animal.effect.GeneSymbol[,c("ID", "F", "p.value","gene.symbol")],
            file="one.way.anova.animal.effect.GeneSymbol.txt", sep="\t", dec=".", row.names=F, quote=F)
#Cleanup
rm(one.way.anova.animal.effect)
rm(one.way.anova.animal.effect.GeneSymbol)


# Plot of top animal-influenced annotated probe set  --------


# Expression values
# Grouping factors (animal + time + treatment)
par(mfrow=c(1,3))
# Plot of expression by animal
boxplot(formula=assayData(farms_informative.65)$exprs["Bt.29815.1.A1_at",]~
          pData(farms_informative.65)$Animal,
        cex.axis=1.5, cex.lab=1.5,
        ylim=c(4, 14),
        xlab="Animal", ylab="log2 (intensity)")

# # Plot of expression by treatment
# boxplot(formula=assayData(farms_informative.65)$exprs["Bt.29815.1.A1_at",]~
#           pData(farms_informative.65)$Treatment, cex.axis=1.5)
# # Plot of expression by timepoint
# boxplot(formula=assayData(farms_informative.65)$exprs["Bt.29815.1.A1_at",]~
#           pData(farms_informative.65)$TimePoint, cex.axis=1.5)
# 
# # Plot of TimePoint/Animal
# boxplot(formula=assayData(farms_informative.65)$exprs["Bt.29815.1.A1_at",]~
#           pData(farms_informative.65)$TimePoint+pData(farms_informative.65)$Animal,
#         las=2, cex.axis=0.8, at=c(1,2,3,4,6,7,8,9,11,12,13,14,16,17,18,19,21,22,23,24))
# # Plot of Treatment/Animal
# boxplot(formula=assayData(farms_informative.65)$exprs["Bt.29815.1.A1_at",]~
#           pData(farms_informative.65)$Treatment+pData(farms_informative.65)$Animal,
#         las=2, cex.axis=0.8, at=c(1,2,3,4,6,7,8,9,11,12,13,14,16,17,18,19,21,22,23,24))
# # Plot of Timepoint/Treatment
# boxplot(formula=assayData(farms_informative.65)$exprs["Bt.29815.1.A1_at",]~
#           pData(farms_informative.65)$Treatment+pData(farms_informative.65)$TimePoint,
#         las=2, cex.axis=0.8, at=c(1,2,3,4,6,7,8,9,11,12,13,14,16,17,18,19))
# Restore the plot space
par(mfrow=c(1,1))

library(ggplot2)

# Bt.4751.1.S1_a_at BOLA-DQA2
d = data.frame(log2Intensity=exprs(farms_informative.65)["Bt.4751.1.S1_a_at",], Animal=pData(farms_informative.65)$Animal)
ggplot(data=d) +
  geom_boxplot(aes(x=Animal, y=log2Intensity, fill=Animal)) +
  ggtitle("Bt.4751.1.S1_a_at = BOLA-DQA2") +
  theme(text = element_text(size=20)) 

# Bt.4751.2.S1_a_at BOLA-DQA2
d = data.frame(log2Intensity=exprs(farms_informative.65)["Bt.4751.2.S1_a_at",], Animal=pData(farms_informative.65)$Animal)
ggplot(data=d) +
  geom_boxplot(aes(x=Animal, y=log2Intensity, fill=Animal)) +
  ggtitle("Bt.4751.2.S1_a_at = BOLA-DQA2") +
  theme(text = element_text(size=20)) 

# Bt.7220.1.S1_at BOLA-DQB
d = data.frame(log2Intensity=exprs(farms_informative.65)["Bt.7220.1.S1_at",], Animal=pData(farms_informative.65)$Animal)
ggplot(data=d) +
  geom_boxplot(aes(x=Animal, y=log2Intensity, fill=Animal)) +
  ggtitle("Bt.7220.1.S1_at = BOLA-DQA2") +
  theme(text = element_text(size=20)) 

# Bt.22867.1.S1_at BOLA-DQA1
d = data.frame(log2Intensity=exprs(farms_informative.65)["Bt.22867.1.S1_at",], Animal=pData(farms_informative.65)$Animal)
ggplot(data=d) +
  geom_boxplot(aes(x=Animal, y=log2Intensity, fill=Animal)) +
  ggtitle("Bt.22867.1.S1_at = BOLA-DQA1") +
  theme(text = element_text(size=20)) 

rm(d)

# Plot of top 40 animal-influenced annotated probe set  --------


# Expression values
par(mfrow=c(1,3))
# Plot of expression by animal
boxplot(formula=assayData(farms_informative.65)$exprs["Bt.17219.1.A1_at",]~
          pData(farms_informative.65)$Animal,
        cex.axis=1.5, cex.lab=1.5,
        ylim=c(4, 14),
        xlab="Animal", ylab="log2 (intensity)")
# Restore the plot space
par(mfrow=c(1,1))


# List of BOVIS specific probe sets at 2 hours ---------------------------------


# Get the list of probe sets DE by BOVIS and MAP at 2 hours
DE.specific.bovis.2hrs = as.character(
  DE.probe.contrast[DE.probe.contrast$"BOVIS_2HR-CN_2HR" == T &
                      DE.probe.contrast$"MAP_2HR-CN_2HR" == F &
                      DE.probe.contrast$"BCG_2HR-CN_2HR" == F,]$ID)
# Write the MAP expression data for those genes
tmp = topTable(ebayes.informative, coef= "BOVIS_2HR-CN_2HR", number=nrow(ebayes.informative), adjust.method="BH")
rownames(tmp) = tmp$ID
write.table(x=tmp[DE.specific.bovis.2hrs,], file="DE.specific.bovis.2hrs_BOVIS.expression.txt", quote=F, sep="\t", row.names=F)
# Cleanup
rm(DE.specific.bovis.2hrs)
rm(tmp)


# List of BOVIS specific probe sets at 6 hours ---------------------------------


# Get the list of probe sets DE by BOVIS and MAP at 6 hours
DE.specific.bovis.6hrs = as.character(
  DE.probe.contrast[DE.probe.contrast$"BOVIS_6HR-CN_6HR" == T &
                      DE.probe.contrast$"MAP_6HR-CN_6HR" == F &
                      DE.probe.contrast$"BCG_6HR-CN_6HR" == F,]$ID)
# Write the MAP expression data for those genes
tmp = topTable(ebayes.informative, coef= "BOVIS_6HR-CN_6HR", number=nrow(ebayes.informative), adjust.method="BH")
rownames(tmp) = tmp$ID
write.table(x=tmp[DE.specific.bovis.6hrs,], file="DE.specific.bovis.6hrs_BOVIS.expression.txt", quote=F, sep="\t", row.names=F)
# Cleanup
rm(DE.specific.bovis.6hrs)
rm(tmp)


# List of BOVIS specific probe sets at 24 hours ---------------------------------


# Get the list of probe sets DE by BOVIS and MAP at 24 hours
DE.specific.bovis.24hrs = as.character(
  DE.probe.contrast[DE.probe.contrast$"BOVIS_24HR-CN_24HR" == T &
                      DE.probe.contrast$"MAP_24HR-CN_24HR" == F &
                      DE.probe.contrast$"BCG_24HR-CN_24HR" == F,]$ID)
# Write the MAP expression data for those genes
tmp = topTable(ebayes.informative, coef= "BOVIS_24HR-CN_24HR", number=nrow(ebayes.informative), adjust.method="BH")
rownames(tmp) = tmp$ID
write.table(x=tmp[DE.specific.bovis.24hrs,], file="DE.specific.bovis.24hrs_BOVIS.expression.txt", quote=F, sep="\t", row.names=F)
# Cleanup
rm(DE.specific.bovis.24hrs)
rm(tmp)


# List of MAP-notBCG probe sets at 2 hours ---------------------------------


# Get the list of probe sets DE by BOVIS and MAP at 2 hours
DE.map.notBCG.2hrs = as.character(
  DE.probe.contrast[DE.probe.contrast$"MAP_2HR-CN_2HR" == T &
                      DE.probe.contrast$"BCG_2HR-CN_2HR" == F,]$ID)
# Write the MAP expression data for those genes
tmp = topTable(ebayes.informative, coef= "MAP_2HR-CN_2HR", number=nrow(ebayes.informative), adjust.method="BH")
rownames(tmp) = tmp$ID
write.table(x=tmp[DE.map.notBCG.2hrs,], file="DE.map.notBCG.2hrs_MAP.expression.txt", quote=F, sep="\t", row.names=F)
# Cleanup
rm(DE.map.notBCG.2hrs)
rm(tmp)


# List of MAP-notBCG probe sets at 6 hours ---------------------------------


# Get the list of probe sets DE by BOVIS and MAP at 6 hours
DE.map.notBCG.6hrs = as.character(
  DE.probe.contrast[DE.probe.contrast$"MAP_6HR-CN_6HR" == T &
                      DE.probe.contrast$"BCG_6HR-CN_6HR" == F,]$ID)
# Write the MAP expression data for those genes
tmp = topTable(ebayes.informative, coef= "MAP_6HR-CN_6HR", number=nrow(ebayes.informative), adjust.method="BH")
rownames(tmp) = tmp$ID
write.table(x=tmp[DE.map.notBCG.6hrs,], file="DE.map.notBCG.6hrs_MAP.expression.txt", quote=F, sep="\t", row.names=F)
# Cleanup
rm(DE.map.notBCG.6hrs)
rm(tmp)


# List of MAP-notBCG probe sets at 24 hours ---------------------------------


# Get the list of probe sets DE by BOVIS and MAP at 24 hours
DE.map.notBCG.24hrs = as.character(
  DE.probe.contrast[DE.probe.contrast$"MAP_24HR-CN_24HR" == T &
                      DE.probe.contrast$"BCG_24HR-CN_24HR" == F,]$ID)
# Write the MAP expression data for those genes
tmp = topTable(ebayes.informative, coef= "MAP_24HR-CN_24HR", number=nrow(ebayes.informative), adjust.method="BH")
rownames(tmp) = tmp$ID
write.table(x=tmp[DE.map.notBCG.24hrs,], file="DE.map.notBCG.24hrs_MAP.expression.txt", quote=F, sep="\t", row.names=F)
# Cleanup
rm(DE.map.notBCG.24hrs)
rm(tmp)


# List of BCG-notMAP probe sets at 2 hours ---------------------------------


# Get the list of probe sets DE by BOVIS and BCG at 2 hours
DE.BCG.notMAP.2hrs = as.character(
  DE.probe.contrast[DE.probe.contrast$"MAP_2HR-CN_2HR" == F &
                      DE.probe.contrast$"BCG_2HR-CN_2HR" == T,]$ID)
# Write the MAP expression data for those genes
tmp = topTable(ebayes.informative, coef= "BCG_2HR-CN_2HR", number=nrow(ebayes.informative), adjust.method="BH")
rownames(tmp) = tmp$ID
write.table(x=tmp[DE.BCG.notMAP.2hrs,], file="DE.BCG.notMAP.2hrs_BCG.expression.txt", quote=F, sep="\t", row.names=F)
# Cleanup
rm(DE.BCG.notMAP.2hrs)
rm(tmp)


# List of BCG-notMAP  probe sets at 6 hours ---------------------------------


# Get the list of probe sets DE by BCG but not MAP at 6 hours
DE.BCG.notMAP.6hrs = as.character(
  DE.probe.contrast[DE.probe.contrast$"MAP_6HR-CN_6HR" == F &
                      DE.probe.contrast$"BCG_6HR-CN_6HR" == T,]$ID)
# Write the MAP expression data for those genes
tmp = topTable(ebayes.informative, coef= "BCG_6HR-CN_6HR", number=nrow(ebayes.informative), adjust.method="BH")
rownames(tmp) = tmp$ID
write.table(x=tmp[DE.BCG.notMAP.6hrs,], file="DE.BCG.notMAP.6hrs_BCG.expression.txt", quote=F, sep="\t", row.names=F)
# Cleanup
rm(DE.BCG.notMAP.6hrs)
rm(tmp)


# List of BCG-notMAP  probe sets at 24 hours ---------------------------------


# Get the list of probe sets DE by BOVIS and BCG at 24 hours
DE.BCG.notMAP.24hrs = as.character(
  DE.probe.contrast[DE.probe.contrast$"MAP_24HR-CN_24HR" == F &
                      DE.probe.contrast$"BCG_24HR-CN_24HR" == T,]$ID)
# Write the MAP expression data for those genes
tmp = topTable(ebayes.informative, coef= "BCG_24HR-CN_24HR", number=nrow(ebayes.informative), adjust.method="BH")
rownames(tmp) = tmp$ID
write.table(x=tmp[DE.BCG.notMAP.24hrs,], file="DE.BCG.notMAP.24hrs_BCG.expression.txt", quote=F, sep="\t", row.names=F)
# Cleanup
rm(DE.BCG.notMAP.24hrs)
rm(tmp)


# List of ENSEMBL-annotated BCG-notMAP  probe sets at 6 hours ---------------------------------


# Get the list of probe sets DE by BCG but not MAP at 6 hours
DE.BCG.notMAP.6hrs = as.character(
  DE.probe.contrast[DE.probe.contrast$"MAP_6HR-CN_6HR" == F &
                      DE.probe.contrast$"BCG_6HR-CN_6HR" == T,]$ID)
# Write the MAP expression data for those genes
tmp = topTable(ebayes.informative, coef= "BCG_6HR-CN_6HR", number=nrow(ebayes.informative), adjust.method="BH")
rownames(tmp) = tmp$ID
write.table(x=AnnotateProbesTable(df=tmp[DE.BCG.notMAP.6hrs,], gene.symbol=F, ENSEMBL=T, ENTREZ=F),
            file="DE.BCG.notMAP.6hrs_ENSEMBL.annotated_BCG.expression.txt", quote=F, sep="\t", row.names=F)
# Cleanup
rm(DE.BCG.notMAP.6hrs)
rm(tmp)


# List of ENSEMBL-annotated BOVIS specific probe sets at 2 hours ---------------------------------

# Get the list of probe sets 
DE.specific.bovis.2hrs = as.character(
  DE.probe.contrast[DE.probe.contrast$"BOVIS_2HR-CN_2HR" == T &
                      DE.probe.contrast$"MAP_2HR-CN_2HR" == F &
                      DE.probe.contrast$"BCG_2HR-CN_2HR" == F,]$ID)
# Write the MAP expression data for those genes
tmp = topTable(ebayes.informative, coef= "BOVIS_2HR-CN_2HR", number=nrow(ebayes.informative), adjust.method="BH")
rownames(tmp) = tmp$ID
write.table(x=AnnotateProbesTable(df=tmp[DE.specific.bovis.2hrs,], gene.symbol=F, ENSEMBL=T, ENTREZ=F),
                                  file="DE.specific.bovis.2hrs_ENSEMBL.annotated_BOVIS.expression.txt", quote=F, sep="\t", row.names=F)
# Cleanup
rm(DE.specific.bovis.2hrs)
rm(tmp)


# List of DE genes shared at 2 hpi (MAP+BCG+BOVI)


# Get the list of probe sets DE common to all at 2 hpi  ---------------------------------
DE.common.2hrs = as.character(
  DE.probe.contrast[DE.probe.contrast$"BOVIS_2HR-CN_2HR" == T &
                      DE.probe.contrast$"MAP_2HR-CN_2HR" == T &
                      DE.probe.contrast$"BCG_2HR-CN_2HR" == T,]$ID)

## The following will gather the data of the three mycobacteria-control contrasts
## in a single dataset. Therefore, cytoscape will be able to color the network nodes
## with either of the contrast results from a single file.

# get the logFC, AveExpr and adj.P.Val for BOVIS-CN contrast 2 hpi
tmp = topTable(ebayes.informative, coef= "BOVIS_2HR-CN_2HR", number=nrow(ebayes.informative), adjust.method="BH")
# Keep only the columns I usually upload in Cytoscape
topTable.merged.2hpi = tmp[,c("ID","logFC","AveExpr","adj.P.Val")]
# Prefix the columns with the specific contrast
colnames(topTable.merged.2hpi) = c("ID", paste("BOV.CN", colnames(topTable.merged.2hpi)[2:4], sep = "."))

# Now get the logFC, AveExpr and adj.P.Val for BCG-CN contrast 2 hpi 
tmp = topTable(ebayes.informative, coef= "BCG_2HR-CN_2HR", number=nrow(ebayes.informative), adjust.method="BH")
rownames(tmp) = tmp$ID
tmp = tmp[,c("ID","logFC","AveExpr","adj.P.Val")]
colnames(tmp) = c("ID", paste("BCG.CN", colnames(tmp)[2:4], sep = "."))
# Merge with the final table
topTable.merged.2hpi = merge(x=topTable.merged.2hpi, y=tmp, by="ID")

# Repeat with the MAP-CN cotnrast at 2 hpi
tmp = topTable(ebayes.informative, coef= "MAP_2HR-CN_2HR", number=nrow(ebayes.informative), adjust.method="BH")
rownames(tmp) = tmp$ID
tmp = tmp[,c("ID","logFC","AveExpr","adj.P.Val")]
colnames(tmp) = c("ID", paste("MAP.CN", colnames(tmp)[2:4], sep = "."))
topTable.merged.2hpi = merge(x=topTable.merged.2hpi, y=tmp, by="ID")

# Restrict to the DE genes common to all mycobacteria
rownames(topTable.merged.2hpi) = topTable.merged.2hpi$ID

# Write the table out with ensembl annotations
write.table(x=AnnotateProbesTable(df=topTable.merged.2hpi[DE.common.2hrs,], gene.symbol=F, ENSEMBL=T, ENTREZ=F),
            file="InnateDB-cytoscape-BINGO-cerebral/common.DE.2hrs_ENSEMBL.annotated_BOVIS.BCG.MAP.expression.txt", quote=F, sep="\t", row.names=F)


head(topTable.merged.2hpi)
# Cleanup
rm(tmp)
rm(DE.common.2hrs)
rm(topTable.merged.2hpi)


# Expression intensities boxplot built-in function ---------------------------------


boxplot_expression= function(eSet, probeset, groupCol = "Group", groups="",
                             annotPkg = "bovine.db", ylim=c(4,16), cex.axis=0.78, las=3, col="", values=FALSE)
{
  # Given an expression set, plots the expression intensities according to a grouping factor.
  #
  # Args:
  #   eSet: expression set
  #   probeset: probe set of which to plot expression
  #   groups: Layout of groups in the plot from left to right. Default is order in pData.
  #   annotPkg:  package where to look for annotation of probeset. Default is bovine.db.
  #   ylim: range of the Y-axis. Default is c(4, 16)
  #   cex.axis: % size of X-axis labels. Default is 0.78.
  #   las:  Orientation of X-axis labels. Default is vertical.
  #   col:  color of boxes. Default is rainbox colormap.
  #   values: Whether to return the values plotted. Default is "FALSE"
  # Returns:
  #     The values plotted, if required.
  #
  # Required packages: Biobase, annotate
  #
  # Packages required to annotate
  library(annotate, quietly=TRUE)
  # Fancy way to pass the library to load as an argument
  library(annotPkg, character.only=TRUE, quietly=TRUE)
  # Error handling
  if (!probeset %in% rownames(assayData(eSet)$exprs)){
    stop("Probeset not found in expression set: ", probeset)
  }
  # Initialises the data frame to plot with condition names
  df = data.frame(row.names=colnames(assayData(eSet)$exprs))
  # Fills the data frame with the expression values
  df[,probeset] = assayData(eSet)$exprs[probeset,]
  # Adds grouping information
  if (any(groups == "")){
    groups = unique(pData(eSet)[,groupCol])
  }
  df$group = factor(pData(eSet)[,groupCol], groups)
  # Colors by group
  n = length(levels(df$group)) # Number of groups = number of colors
  print(levels(df$group))
  if (any(col != "") & length(col) != n){
    stop("Number of colors is different from number of groups: colors(", 
         length(col), ") and groups(", n,").")
  }
  if (any(col == "")){
    col = rainbow(n, s = 1, v = 1, start = 0, end = max(1, n - 1)/n, alpha = 1)
  }
  # Obtains the gene symbol annotations
  gene.symbol = getSYMBOL(probeset, annotPkg)
  # plots the data
  boxplot(df[,probeset]~group, data=df, ylim = ylim,
          cex.axis=cex.axis, las=las,
          main=paste(probeset, "=", gene.symbol, "| Expression plot"),
          col=col)
  # If requested return the values plotted
  if(values){
    return(df)
  }
}


# Expression boxplots using built-in function --------------------------------------

library(Cairo)
# IL1B
Cairo(width=800, height=800, file="Expression.boxplots/IL1B.Bt.4856.1.S1_at.png", type="png", bg="white")
boxplot_expression(eSet=farms_informative.65,
                   probeset="Bt.4856.1.S1_at",
                   groupCol="Group",
                   groups=c("CN_0HR","CN_2HR","CN_6HR","CN_24HR",
                            "MAP_2HR","MAP_6HR","MAP_24HR",
                            "BCG_2HR","BCG_6HR","BCG_24HR",
                            "BOVIS_2HR","BOVIS_6HR","BOVIS_24HR"),
                   annotPkg="bovine.db",
                   ylim=c(4,14),
                   cex.axis=0.7,
                   las=2,
                   col=c(rep("lightblue",4),rep("lightgreen",3),rep("gold",3),rep("plum",3)),
                   values=F)
dev.off()

# IFNAC
Cairo(width=800, height=800, file="Expression.boxplots/IFNAC.Bt.416.1.S1_at.png", type="png", bg="white")
boxplot_expression(eSet=farms_informative.65,
                   probeset="Bt.416.1.S1_at",
                   groupCol="Group",
                   groups=c("CN_0HR","CN_2HR","CN_6HR","CN_24HR",
                            "MAP_2HR","MAP_6HR","MAP_24HR",
                            "BCG_2HR","BCG_6HR","BCG_24HR",
                            "BOVIS_2HR","BOVIS_6HR","BOVIS_24HR"),
                   annotPkg="bovine.db",
                   ylim=c(6,7),
                   cex.axis=0.7,
                   las=2,
                   col=c(rep("lightblue",4),rep("lightgreen",3),rep("gold",3),rep("plum",3)),
                   values=F)
dev.off()

# IFNG
Cairo(width=800, height=800, file="Expression.boxplots/IFNG.Bt.188.1.S1_at.png", type="png", bg="white")
boxplot_expression(eSet=farms_informative.65,
                   probeset="Bt.188.1.S1_at",
                   groupCol="Group",
                   groups=c("CN_0HR","CN_2HR","CN_6HR","CN_24HR",
                            "MAP_2HR","MAP_6HR","MAP_24HR",
                            "BCG_2HR","BCG_6HR","BCG_24HR",
                            "BOVIS_2HR","BOVIS_6HR","BOVIS_24HR"),
                   annotPkg="bovine.db",
                   ylim=c(4,14),
                   cex.axis=0.7,
                   las=2,
                   col=c(rep("lightblue",4),rep("lightgreen",3),rep("gold",3),rep("plum",3)),
                   values=F)
dev.off()


# Expression intensities boxplot ggplot2 function (Doesn't work) -------------------------

gg.boxplot_expression= function(eSet, probeset, groupCol = "Group", groups="",
                             annotPkg = "bovine.db", ylim=c(4,16), cex.axis=0.78, las=3, col="", values=FALSE)
{
  # Given an expression set, plots the expression intensities according to a grouping factor.
  #
  # Args:
  #   eSet: expression set
  #   probeset: probe set of which to plot expression
  #   groups: Layout of groups in the plot from left to right. Default is order in pData.
  #   annotPkg:  package where to look for annotation of probeset. Default is bovine.db.
  #   ylim: range of the Y-axis. Default is c(4, 16)
  #   cex.axis: % size of X-axis labels. Default is 0.78.
  #   las:  Orientation of X-axis labels. Default is vertical.
  #   col:  color of boxes. Default is rainbox colormap.
  #   values: Whether to return the values plotted. Default is "FALSE"
  # Returns:
  #     The values plotted, if required.
  #
  # Required packages: Biobase, annotate
  #
  # Packages required to annotate
  library(annotate, quietly=TRUE)
  # Fancy way to pass the library to load as an argument
  library(annotPkg, character.only=TRUE, quietly=TRUE)
  # Plotting library
  library(ggplot2, quietly=TRUE)
  # Error handling
  if (!probeset %in% rownames(assayData(eSet)$exprs)){
    stop("Probeset not found in expression set: ", probeset)
  }
  # Initialises the data frame to plot with condition names
  df = data.frame(row.names=colnames(assayData(eSet)$exprs))
  # Fills the data frame with the expression values
  df[,probeset] = assayData(eSet)$exprs[probeset,]
  # Adds grouping information
  if (any(groups == "")){
    groups = unique(pData(eSet)[,groupCol])
  }
  df$group = factor(pData(eSet)[,groupCol], groups)
  # Colors by group
  n = length(levels(df$group)) # Number of groups = number of colors
  if (any(col != "") & length(col) != n){
    stop("Number of colors is different from number of groups: colors(", 
         length(col), ") and groups(", n,").")
  }
  if (any(col == "")){
    col = rainbow(n, s = 1, v = 1, start = 0, end = max(1, n - 1)/n, alpha = 1)
  }
  # Obtains the gene symbol annotations
  gene.symbol = getSYMBOL(probeset, annotPkg)
  # Data-driven ggplot definition
  tmp.graph = ggplot(df, aes(x=df$group,
                             y=df[,probeset]))
  # Layout-driven ggplot definition
  tmp.graph + geom_boxplot(fill=col) +
    ggtitle(paste(probeset, "=", gene.symbol, "| Expression plot")) +
    theme(text = element_text(size=20), 
          axis.text.x = element_text(angle=90, face="bold"),
          plot.title = element_text(vjust=1.5),
          axis.title.x = element_text(vjust=-0.3),
          axis.title.y = element_text(vjust=0.25)
    ) + 
    xlab("") +
    ylab("log2 intensity") + 
    ylim(ylim[1], ylim[2])
  # If requested return the values plotted
  if(values){
    return(df)
  }
}


# Expression boxplots using ggplot2 function ------------------------------

library(Cairo)

##########
## IL1B ##
##########
# define the probeset to plot
tmp.probeset = "Bt.4856.1.S1_at" 
# creates a dataset with rows named
tmp = data.frame(row.names=colnames(assayData(farms_informative.65)$exprs))
# creates a column named after the probe ID and fills in with expression data
tmp[,tmp.probeset] = assayData(farms_informative.65)$exprs[tmp.probeset,]
# creates a column with grouping factor based on treatment and time
tmp$group = factor(pData(farms_informative.65)[,"Group"], c("CN_0HR","CN_2HR","CN_6HR","CN_24HR","MAP_2HR","MAP_6HR","MAP_24HR","BCG_2HR","BCG_6HR","BCG_24HR","BOVIS_2HR","BOVIS_6HR","BOVIS_24HR"))
# Prepares the ggplot
tmp.graph = ggplot(data=tmp, aes(x=tmp$group, y=tmp[,tmp.probeset]))
# Opens a Cairo quality file
Cairo(width=800, height=800, file="Expression.boxplots/IL1B.Bt.4856.1.S1_at.(ggplot2).png", type="png", bg="white")
# Outputs the graph in Cairo
tmp.graph +
  ggtitle(tmp.probeset) +
  geom_boxplot(fill=c(rep("lightblue",4),rep("lightgreen",3),rep("gold",3),rep("plum",3))) +
  theme(text = element_text(size=20), 
        axis.text.x = element_text(vjust=0.5, angle=90, face="bold"),
        plot.title = element_text(vjust=1.5),
        axis.title.x = element_text(vjust=-0.3), # no axis title anyway
        axis.title.y = element_text(vjust=0.25)
  ) + 
  xlab("") +
  ylab("") + 
  ylim(4, 14)
# Writes the file to disk
dev.off()
#Cleanup
rm(tmp, tmp.graph, tmp.probeset)

###########
## IFNAC ##
###########
# define the probeset to plot
tmp.probeset = "Bt.416.1.S1_at"
# creates a dataset with rows named
tmp = data.frame(row.names=colnames(assayData(farms_informative.65)$exprs))
# creates a column named after the probe ID and fills in with expression data
tmp[,tmp.probeset] = assayData(farms_informative.65)$exprs[tmp.probeset,]
# creates a column with grouping factor based on treatment and time
tmp$group = factor(pData(farms_informative.65)[,"Group"], c("CN_0HR","CN_2HR","CN_6HR","CN_24HR","MAP_2HR","MAP_6HR","MAP_24HR","BCG_2HR","BCG_6HR","BCG_24HR","BOVIS_2HR","BOVIS_6HR","BOVIS_24HR"))
# Prepares the ggplot
tmp.graph = ggplot(data=tmp, aes(x=tmp$group, y=tmp[,tmp.probeset]))
# Opens a Cairo quality file
Cairo(width=800, height=800, file="Expression.boxplots/IFNAC.Bt.416.1.S1_at.(ggplot2).png", type="png", bg="white")
# Outputs the graph in Cairo
tmp.graph +
  ggtitle(tmp.probeset) +
  geom_boxplot(fill=c(rep("lightblue",4),rep("lightgreen",3),rep("gold",3),rep("plum",3))) +
  theme(text = element_text(size=20), 
        axis.text.x = element_text(vjust=0.5, angle=90, face="bold"),
        plot.title = element_text(vjust=1.5),
        axis.title.x = element_text(vjust=-0.3), # no axis title anyway
        axis.title.y = element_text(vjust=0.25)
  ) + 
  xlab("") +
  ylab("") + 
  ylim(6, 7)
# Writes the file to disk
dev.off()
#Cleanup
rm(tmp, tmp.graph, tmp.probeset)

##########
## IFNG ##
##########
# define the probeset to plot
tmp.probeset = "Bt.188.1.S1_at"
# creates a dataset with rows named
tmp = data.frame(row.names=colnames(assayData(farms_informative.65)$exprs))
# creates a column named after the probe ID and fills in with expression data
tmp[,tmp.probeset] = assayData(farms_informative.65)$exprs[tmp.probeset,]
# creates a column with grouping factor based on treatment and time
tmp$group = factor(pData(farms_informative.65)[,"Group"], c("CN_0HR","CN_2HR","CN_6HR","CN_24HR","MAP_2HR","MAP_6HR","MAP_24HR","BCG_2HR","BCG_6HR","BCG_24HR","BOVIS_2HR","BOVIS_6HR","BOVIS_24HR"))
# Prepares the ggplot
tmp.graph = ggplot(data=tmp, aes(x=tmp$group, y=tmp[,tmp.probeset]))
# Opens a Cairo quality file
Cairo(width=800, height=800, file="Expression.boxplots/IFNG.Bt.188.1.S1_at.(ggplot2).png", type="png", bg="white")
# Outputs the graph in Cairo
tmp.graph +
  ggtitle(tmp.probeset) +
  geom_boxplot(fill=c(rep("lightblue",4),rep("lightgreen",3),rep("gold",3),rep("plum",3))) +
  theme(text = element_text(size=20), 
        axis.text.x = element_text(vjust=0.5, angle=90, face="bold"),
        plot.title = element_text(vjust=1.5),
        axis.title.x = element_text(vjust=-0.3), # no axis title anyway
        axis.title.y = element_text(vjust=0.25)
  ) + 
  xlab("") +
  ylab("") + 
  ylim(4, 14)
# Writes the file to disk
dev.off()
#Cleanup
rm(tmp, tmp.graph, tmp.probeset)


# Useful variable in this analysis ----------------------------------------


# each column a sample, each row a probe / use t() to transform
head(assayData(farms_informative.65)$exprs)
# sample names (same order as pData)
rownames(t(assayData(farms_informative.65)$exprs))
# Animal for each sample
pData(farms_informative.65)$Animal
# TimePoint for each sample
pData(farms_informative.65)$TimePoint
# Treatment for each sample
pData(farms_informative.65)$Treatment


# Boolean data frame of whether probes up or down DE in contrasts --------


# Function returning whether each probe is DE in each contrast
OppositeDE <- function(ebayes, contrasts, probesets)
{
  DE.table = data.frame(ID=probesets)
  for(contrast in contrasts){
    topT = topTable(ebayes, coef= contrast, number=nrow(ebayes.informative), adjust.method="BH")
    # Checks whether each genes is DE and up or down for that contrast
    topT$bool = "NotDE"
    topT$bool[topT$adj.P.Val < 0.05 & topT$logFC > 0] = "Up"
    topT$bool[topT$adj.P.Val < 0.05 & topT$logFC < 0] = "Down"
    # Renames the column with the contrast name
    names(topT)[names(topT)=="bool"] = contrast
    # put data in data frame
    DE.table = merge(x=DE.table, y=topT[c("ID",contrast)], by="ID", sort=TRUE)
  }
  return(DE.table)
}
# Return whether each probe is DE in each contrast
DE.direction.probe.contrast = OppositeDE(ebayes=ebayes.informative,contrasts=contrasts.DGE.names,
                                probesets=rownames(assayData(farms_informative.65)$exprs))
save(DE.direction.probe.contrast, file="DE.direction.probe.contrast.RData")
write.table(x=DE.direction.probe.contrast, file="DE.direction.probe.contrast.txt", quote=F, sep="\t", row.names=F)

# Cleanup
rm(DE.direction.probe.contrast)


# Look for genes DE up by one TB and down by another ----------------------

load("DE.direction.probe.contrast.RData")

# Subset for the myco-CN at 2 hpi plus the ID column
myco.CN.2HR = DE.direction.probe.contrast[,c(1,grep(pattern="CN_2HR", x=colnames(DE.direction.probe.contrast)))]

# Which probes have Up and Down by different mycobacteria?
# If all contrasts are Up or NotDE, consensus is Up
# If all contrasts are Down or NotDE, consensus is Down
# Otherwise, consensus is "Diff"
myco.CN.2HR$consensus = "Diff"
myco.CN.2HR$consensus[apply(X=myco.CN.2HR[,2:4], MARGIN=1, FUN=function(x) {all(x == "NotDE")})] = "NotDE"
myco.CN.2HR$consensus[apply(X=myco.CN.2HR[,2:4], MARGIN=1, FUN=function(x) {all(x %in% c("Up","NotDE")) & any(x == "Up")})] = "Up"
myco.CN.2HR$consensus[apply(X=myco.CN.2HR[,2:4], MARGIN=1, FUN=function(x) {all(x %in% c("Down","NotDE")) & any(x == "Down")})] = "Down"
# Number of Rows where consensus is Diff ? 0
sum(myco.CN.2HR$consensus == "Diff")
# No opposite trends at 2 hpi

#Cleanup
rm(myco.CN.2HR)


# Subset for the myco-CN at 6 hpi plus the ID column
myco.CN.6HR = DE.direction.probe.contrast[,c(1,grep(pattern="CN_6HR", x=colnames(DE.direction.probe.contrast)))]

# Which probes have Up and Down by different mycobacteria?
# If all contrasts are Up or NotDE, consensus is Up
# If all contrasts are Down or NotDE, consensus is Down
# Otherwise, consensus is "Diff"
myco.CN.6HR$consensus = "Diff"
myco.CN.6HR$consensus[apply(X=myco.CN.6HR[,2:4], MARGIN=1, FUN=function(x) {all(x == "NotDE")})] = "NotDE"
myco.CN.6HR$consensus[apply(X=myco.CN.6HR[,2:4], MARGIN=1, FUN=function(x) {all(x %in% c("Up","NotDE")) & any(x == "Up")})] = "Up"
myco.CN.6HR$consensus[apply(X=myco.CN.6HR[,2:4], MARGIN=1, FUN=function(x) {all(x %in% c("Down","NotDE")) & any(x == "Down")})] = "Down"
# Number of Rows where consensus is Diff ? 0
sum(myco.CN.6HR$consensus == "Diff")
# No opposite trends at 6 hpi

#Cleanup
rm(myco.CN.6HR)


# Subset for the myco-CN at 24 hpi plus the ID column
myco.CN.24HR = DE.direction.probe.contrast[,c(1,grep(pattern="CN_24HR", x=colnames(DE.direction.probe.contrast)))]

# Which probes have Up and Down by different mycobacteria?
# If all contrasts are Up or NotDE, consensus is Up
# If all contrasts are Down or NotDE, consensus is Down
# Otherwise, consensus is "Diff"
myco.CN.24HR$consensus = "Diff"
myco.CN.24HR$consensus[apply(X=myco.CN.24HR[,2:4], MARGIN=1, FUN=function(x) {all(x == "NotDE")})] = "NotDE"
myco.CN.24HR$consensus[apply(X=myco.CN.24HR[,2:4], MARGIN=1, FUN=function(x) {all(x %in% c("Up","NotDE")) & any(x == "Up")})] = "Up"
myco.CN.24HR$consensus[apply(X=myco.CN.24HR[,2:4], MARGIN=1, FUN=function(x) {all(x %in% c("Down","NotDE")) & any(x == "Down")})] = "Down"
# Number of Rows where consensus is Diff ? 0
sum(myco.CN.24HR$consensus == "Diff")
# No opposite trends at 24 hpi

#Cleanup
rm(myco.CN.24HR)


# MDS using informative probesets -----------------------------------------

# Calculate the MDS
MDS_no0H = plotMDS(farms_informative.65_no0H, top=nrow(farms_informative.65_no0H))
# Check that the order of samples is respected
all(names(MDS_no0H$x) == rownames(pData(targets.65_no0H)))
varLabels(targets.65_no0H)
# Prepare the color and character coding of plotted points and legend
Legend = data.frame(labels=unique(pData(targets.65_no0H)$Group))
# 12 colors to use for the different TimePoint-Treatment groups
Legend$col = c("grey80", "grey40", "black",
               "#FFAAAA", "#FF5555", "#FF0000",
               "olivedrab1", "olivedrab3", "olivedrab4",
               "#BDC9EE","#7B94DD","#3A5FCD")
Legend$pch = c(rep(c(17,15,18),4))
Legend$cex = c(rep(c(1.5,1.5,2),4))
Plot = data.frame(sample=rownames(pData(targets.65_no0H)))
Plot$pch = c(17,15,18)[as.numeric(factor(pData(targets.65_no0H)$TimePoint, levels=c("2HR","6HR","24HR"), ordered=T))]
Plot$cex = c(1.5,1.5,2)[as.numeric(factor(pData(targets.65_no0H)$TimePoint, levels=c("2HR","6HR","24HR"), ordered=T))]
Plot$col = col.legend[as.numeric(factor(x=pData(targets.65_no0H)$Group,
                                        levels=unique(pData(targets.65_no0H)$Group),
                                        ordered=T))]


plot(MDS_no0H$x, MDS_no0H$y, col=Plot$col, pch=Plot$pch, cex=Plot$cex, #xlim=c(-1, 0.8), ylim=c(-0.9, 0.9),
     main="MDS informative", xlab="Dimension 1", ylab="Dimension 2")
legend(x="topleft",col=Legend$col, pch=Legend$pch, inset = 0.01, cex=0.90,
       legend=Legend$labels)

# Cleanup
rm(Legend, Plot, MDS_no0H)


# MDS using informative probesets without control samples -------------------------

farms_informative_no_ctrl = farms_informative.65[,-grep(pattern="CN", x=rownames(pData(targets.65)))]
targets_no_ctrl = targets.65[-grep(pattern="CN", x=rownames(pData(targets.65))),]


# Calculate the MDS
MDS_no_ctrl = plotMDS(farms_informative_no_ctrl, top=nrow(farms_informative_no_ctrl))
# Check that the order of samples is respected
all(names(MDS_no0H$x) == rownames(pData(targets.65_no0H)))
varLabels(targets.65_no0H)
# Prepare the color and character coding of plotted points and legend
Legend = data.frame(labels=unique(pData(targets_no_ctrl)$Group))
# 9 colors to use for the different TimePoint-Treatment groups
Legend$col = c("#FFAAAA", "#FF5555", "#FF0000",
               "olivedrab1", "olivedrab3", "olivedrab4",
               "#BDC9EE","#7B94DD","#3A5FCD")
Legend$pch = c(rep(c(17,15,18),3))
Legend$cex = c(rep(c(1.5,1.5,2),3))
Plot = data.frame(sample=rownames(pData(targets_no_ctrl)))
Plot$pch = c(17,15,18)[as.numeric(factor(pData(targets_no_ctrl)$TimePoint, levels=c("2HR","6HR","24HR"), ordered=T))]
Plot$cex = c(1.5,1.5,2)[as.numeric(factor(pData(targets_no_ctrl)$TimePoint, levels=c("2HR","6HR","24HR"), ordered=T))]
Plot$col = Legend$col[as.numeric(factor(x=pData(targets_no_ctrl)$Group,
                                        levels=unique(pData(targets_no_ctrl)$Group),
                                        ordered=T))]


plot(MDS_no_ctrl$x, MDS_no_ctrl$y, col=Plot$col, pch=Plot$pch, cex=Plot$cex, #xlim=c(-1, 0.8), ylim=c(-0.9, 0.9),
     main="MDS informative", xlab="Dimension 1", ylab="Dimension 2")
legend(x="topleft",col=Legend$col, pch=Legend$pch, inset = 0.01, cex=0.90,
       legend=Legend$labels)

# Cleanup
rm(Legend, Plot, MDS_no0H)


# Semi-automated MDS plot using Manhattan distance --------------------------------

# The dist function calculates distance between the rows of a matrix.
# we need to transform the data.frame and save it in a matrix
m = matrix(t(exprs(farms_informative_no_ctrl)),
           nrow=ncol(exprs(farms_informative_no_ctrl)),
           ncol=nrow(exprs(farms_informative_no_ctrl)),
           dimnames=list(colnames(farms_informative_no_ctrl),
                         rownames(farms_informative_no_ctrl)))
#
d = dist(x=m, method="manhattan")
f = cmdscale(d, eig=TRUE, k=2) # k is the number of dim

all(rownames(f$points) == rownames(pData(targets_no_ctrl)))
# Prepare the color and character coding of plotted points and legend
Legend = data.frame(labels=unique(pData(targets_no_ctrl)$Group))
# 9 colors to use for the different TimePoint-Treatment groups
Legend$col = c("#BDC9EE","#7B94DD","#3A5FCD",
               "olivedrab1", "olivedrab3", "olivedrab4",
               "#FFAAAA", "#FF5555", "#FF0000")
Legend$pch = c(rep(c(17,15,18),3))
Legend$cex = c(rep(c(1.5,1.5,2),3))
Plot = data.frame(sample=rownames(pData(targets_no_ctrl)))
Plot$pch = c(17,15,18)[as.numeric(factor(pData(targets_no_ctrl)$TimePoint, levels=c("2HR","6HR","24HR"), ordered=T))]
Plot$cex = c(1.5,1.5,2)[as.numeric(factor(pData(targets_no_ctrl)$TimePoint, levels=c("2HR","6HR","24HR"), ordered=T))]
Plot$col = Legend$col[as.numeric(factor(x=pData(targets_no_ctrl)$Group,
                                        levels=unique(pData(targets_no_ctrl)$Group),
                                        ordered=T))]

plot(f$points[,1], f$points[,2], xlab="Coordinate 1", ylab="Coordinate 2", 
     col=Plot$col, pch=Plot$pch, cex=Plot$cex)
legend(x="topright",col=Legend$col, pch=Legend$pch, inset = 0.01, cex=0.90,
       legend=Legend$labels)
identify(f$points[,1], f$points[,2], labels = sapply(strsplit(row.names(m), "_"), "[", 1), cex=.7)
#plot(f$points[,1], f$points[,2], xlab="Coordinate 1", ylab="Coordinate 2", 
#     main="Metric  MDS",	 type="n")
#

#Cleanup
rm(m, d, f, Legend, Plot)


# Semi-automated MDS plot using Manhattan distance and ggplot2 ---------------

# The dist function calculates distance between the rows of a matrix.
# we need to transform the data.frame and save it in a matrix
m = matrix(t(exprs(farms_informative_no_ctrl)),
           nrow=ncol(exprs(farms_informative_no_ctrl)),
           ncol=nrow(exprs(farms_informative_no_ctrl)),
           dimnames=list(colnames(farms_informative_no_ctrl),
                         rownames(farms_informative_no_ctrl)))
#
d = dist(x=m, method="manhattan")
f = cmdscale(d, eig=TRUE, k=2) # k is the number of dim

# Check
all(rownames(f$points) == rownames(pData(targets_no_ctrl)))

library(ggplot2)

mds = data.frame(f$points)
# Coordinates are named X1 and X2

rownames(mds)
# we can identify the animal, treatment and time from the rowname
mds$Animal = sapply(X=rownames(mds), FUN=function(x){strsplit(x, "_")[[1]][[1]]})
mds$Treatment = sapply(X=rownames(mds), FUN=function(x){strsplit(x, "_")[[1]][[2]]})
mds$Time = sapply(X=rownames(mds), FUN=function(x){strsplit(x, "_")[[1]][[3]]})
mds$Time[mds$Time == "25HR"] = "24HR"

mds$Time = factor(mds$Time, levels=c("2HR","6HR","24HR"), ordered=TRUE)

ggplot(data=mds) +
  geom_point(aes(x=X1, y=X2,
                 colour=Treatment,
                 shape=Time),
             size = 5) +
  xlab("Dimension 1") +
  ylab("Dimension 2") +
  geom_text(aes(x=X1, y=X2, label=Animal), hjust=1, vjust=-1)



#Cleanup
rm(m, d, f)


# Generate the table of log2FC for given probesets ------------------------

generate_tableFC = function(probesets){
  # Depends on the variable FC.probe.contrast being in the workspace
  # Depends on the variable FC.probe.contrast being in the workspace
  # Annotate the fold change table with gene symbols
  annotated_FC = AnnotateProbesTable(df=FC.probe.contrast, ENSEMBL=F, ENTREZ=F, probeCol="ID")
  # Reorder and subset the fold change table to the desired format to minimise work on Excel later
  annotated_FC= annotated_FC[,c("gene.symbol","BOVIS_2HR-CN_2HR","BOVIS_6HR-CN_6HR","BOVIS_24HR-CN_24HR",
                                "BCG_2HR-CN_2HR","BCG_6HR-CN_6HR","BCG_24HR-CN_24HR",
                                "MAP_2HR-CN_2HR","MAP_6HR-CN_6HR","MAP_24HR-CN_24HR",
                                "BOVIS_2HR-BCG_2HR","BOVIS_6HR-BCG_6HR","BOVIS_24HR-BCG_24HR",
                                "BOVIS_2HR-MAP_2HR","BOVIS_6HR-MAP_6HR","BOVIS_24HR-MAP_24HR")]
  # Extract the rows of interest
  extract_FC = annotated_FC[probesets,]
  # Reorder and subset the DE boolean table to the desired format to minimise work on Excel later
  tmpDE = DE.probe.contrast[,c("ID", "BOVIS_2HR-CN_2HR","BOVIS_6HR-CN_6HR","BOVIS_24HR-CN_24HR",
                                "BCG_2HR-CN_2HR","BCG_6HR-CN_6HR","BCG_24HR-CN_24HR",
                                "MAP_2HR-CN_2HR","MAP_6HR-CN_6HR","MAP_24HR-CN_24HR",
                                "BOVIS_2HR-BCG_2HR","BOVIS_6HR-BCG_6HR","BOVIS_24HR-BCG_24HR",
                                "BOVIS_2HR-MAP_2HR","BOVIS_6HR-MAP_6HR","BOVIS_24HR-MAP_24HR")]
  # Extract the rows of interest
  extract_DE = tmpDE[probesets,]
  # Hack to avoid replacing the gene.symbol by 0
  extract_DE$ID = TRUE
  #
  replacedFC = replace(extract_FC, extract_DE == 0, 0)
  return(replacedFC)
}

curated_table_FC = generate_tableFC(c("Bt.4856.1.S1_at","Bt.191.1.S1_at", # IL-1 signalling
                   "BtAffx.1.7.S1_at","Bt.416.1.S1_at","Bt.16966.1.S1_at","Bt.24795.1.A1_at","Bt.4675.1.S1_a_at","Bt.8143.1.S1_at","Bt.4557.1.S1_at", "Bt.5508.1.S1_at", # IFN-I
                   "Bt.188.1.S1_at", "Bt.4251.1.S1_at", "Bt.4251.2.A1_at", # IFN-II
                   "Bt.26983.1.S1_at", # CASPASES
                   "Bt.9309.1.A1_at", "Bt.20288.1.S1_at", "Bt.12756.1.S1_at", "Bt.3686.1.S1_at", # NFkB
                   "Bt.2899.1.S1_at" # dow-regulated (FOS)
                   ))
write.table(x=curated_table_FC, file="Heatmap.curated.genes.table/curated_table_FC.txt", quote=F, sep="\t", row.names=T, col.names=T)


# Save the expression data per sample -------------------------------------
exprs_farms_informative65 = exprs(farms_informative.65)
save(exprs_farms_informative65, file="exprs_farms_informative65.RData")
save(targets.65, file="targets65.RData")


# Review figure: probeset with variance killing DE calls ------------------

# After-note: remember that DE genes between MAP and BCG are only found at 24 hpi
# and those are only 4 genes.

# Identify probesets DE in BCG-CN but not BCG-MAP nor MAP-CN
BCG.CN.spe = as.character(
  DE.probe.contrast[DE.probe.contrast$"MAP_6HR-CN_6HR" == F &
                      DE.probe.contrast$"BCG_6HR-CN_6HR" == T &
                      DE.probe.contrast$"MAP_6HR-BCG_6HR" == F,]$ID)

# Create a dataset without BOVIS treatment but with CN
farms_informative_no_bovis = farms_informative.65_no0H[,targets.65_no0H$Treatment != "BOVIS"]
targets_no_bovis = targets.65_no0H[targets.65_no0H$Treatment != "BOVIS",]

# 
library(biomaRt)
mart = useMart("ensembl", "btaurus_gene_ensembl")

setwd("C:/Users/krue/Documents/Kevin-Logs/20130320_MDM/Writing thesis chapter")

# Sample one probeset randomly
tmp.probeset = sample(BCG.CN.spe, 1)

# Extract the expression level of that probeset in CN, MAP and BCG in samples from 6 hpi
d = data.frame(log2Intensity=exprs(farms_informative_no_bovis)[tmp.probeset,targets_no_bovis$TimePoint == "6HR"])

# From the rowname we can guess which treatment each sample is from
d$Treatment = sapply(X=rownames(d), FUN=function(x){strsplit(x, "_")[[1]][[2]]})
d$Treatment = factor(d$Treatment, levels=c("CN", "MAP", "BCG"), ordered=TRUE)

gg = ggplot(data=d) +
  geom_boxplot(aes(x=Treatment, y=log2Intensity,
                   fill=Treatment)) +
  ggtitle(paste(tmp.probeset, "=", getBM(attributes=c("external_gene_id"), filters="affy_bovine", values=tmp.probeset, mart=mart)[[1]])) +
  ylab("log2(intensity)") +
  theme(text = element_text(size=20),
        axis.text.x = element_text(size=15, vjust=1)) 

ggsave(filename=paste(tmp.probeset, ".pdf", sep=""), plot=gg)

rm(gg, d, tmp.probeset, BCG.CN.spe)


# Same process except now we search for a gene which is DE in all three contrasts

# Identify probesets DE in BCG-CN, BCG-MAP and MAP-CN (none)
MAP_BCG.shared = as.character(
  DE.probe.contrast[DE.probe.contrast$"MAP_6HR-CN_6HR" == T &
                      DE.probe.contrast$"BCG_6HR-CN_6HR" == T &
                      DE.probe.contrast$"MAP_6HR-BCG_6HR" == F,]$ID)

# Sample one probeset randomly
tmp.probeset = sample(MAP_BCG.shared, 1)

# Extract the expression level of that probeset in CN, MAP and BCG in samples from 6 hpi
d = data.frame(log2Intensity=exprs(farms_informative_no_bovis)[tmp.probeset,targets_no_bovis$TimePoint == "6HR"])

# From the rowname we can guess which treatment each sample is from
d$Treatment = sapply(X=rownames(d), FUN=function(x){strsplit(x, "_")[[1]][[2]]})
d$Treatment = factor(d$Treatment, levels=c("CN", "MAP", "BCG"), ordered=TRUE)

gg = ggplot(data=d) +
  geom_boxplot(aes(x=Treatment, y=log2Intensity,
                   fill=Treatment)) +
  ggtitle(paste(tmp.probeset, "=", getBM(attributes=c("external_gene_id"), filters="affy_bovine", values=tmp.probeset, mart=mart)[[1]])) +
  ylab("log2(intensity)") +
  theme(text = element_text(size=20),
        axis.text.x = element_text(size=15, vjust=1)) 

ggsave(filename=paste(tmp.probeset, "_MAPandBCG.pdf", sep=""), plot=gg)

rm(gg, d, tmp.probeset, MAP_BCG.shared)

# Save expression data for Dave to do his own boxplots
d = data.frame(log2Intensity=exprs(farms_informative.65_no0H)["Bt.24741.1.S1_at",targets.65_no0H$TimePoint == "6HR"])
write.table(x=d, file="C:/Users/krue/Documents/Kevin-Logs/20130320_MDM/Writing thesis chapter/Expression.ggplot/Bt.24741.1.S1_at_no0H.txt", quote=F, sep="\t")
d = data.frame(log2Intensity=exprs(farms_informative.65_no0H)["Bt.5420.1.A1_at",targets.65_no0H$TimePoint == "6HR"])
write.table(x=d, file="C:/Users/krue/Documents/Kevin-Logs/20130320_MDM/Writing thesis chapter/Expression.ggplot/Bt.5420.1.A1_at-no0H.txt", quote=F, sep="\t")

# Dave wants the p-values too to put stars in the significant contrasts


library(limma) # topTable()
# Function returning the fold change in each contrast
extractPValues <- function(ebayes, contrasts, probesets)
{
  p.table = data.frame(ID=probesets)
  for(contrast in contrasts){
    topT = topTable(ebayes, coef=contrast, number=nrow(ebayes), adjust.method="BH")
    # Renames the column with the contrast name
    names(topT)[names(topT)=="adj.P.Val"] = contrast
    topT$ID = rownames(topT)
    # put data in data frame
    p.table = merge(x=p.table, y=topT[c("ID", contrast)], by="ID", sort=TRUE)
  }
  return(p.table)
}
# Obtains the fold change values for each contrast
FDR.probe.contrast = extractPValues(ebayes=ebayes.informative, contrasts=contrasts.DGE.names,
                                       probesets=rownames(assayData(farms_informative.65)$exprs))

save(FDR.probe.contrast, file="FDR.probe.contrast.RData")

# Cleanup
rm(FDR.probe.contrast)


# The p values of the probesets that Dave actually wants
d = extractPValues(ebayes=ebayes.informative, contrasts=contrasts.DGE.names,
               probesets=c("Bt.24741.1.S1_at","Bt.5420.1.A1_at"))
rownames(d) = d$ID
d$ID = NULL



write.table(t(d), file="Expression.ggplot/adj.P.Values.txt", quote=F)


# Export the topTable for each contrast -----------------------------------


library(limma) # topTable()

topTable.dir = "topTables"
dir.create(topTable.dir)

for(contrast in contrasts.DGE.names){
  topT = topTable(ebayes.informative, coef= contrast, number=nrow(ebayes.informative), adjust.method="BH")
  write.table(topT, file = file.path(topTable.dir, paste(contrast, "topTable.txt", sep = ".")), quote = F, sep = "\t", row.names = T, col.names = T)
  save(topT, file = file.path(topTable.dir, paste(contrast, "topTable.RData", sep = ".")))
}
rm(topT)

for(contrast in contrasts.DGE.names){
  assign(x = contrast, value = topTable(ebayes.informative, coef= contrast, number=nrow(ebayes.informative), adjust.method="BH"))
  write.table(get(contrast), file = file.path(topTable.dir, paste(contrast, "topTable.txt", sep = ".")), quote = F, sep = "\t", row.names = T, col.names = T)
  save(list = contrast, file = file.path(topTable.dir, paste(contrast, "topTable.RData", sep = ".")))
}
rm(contrast)

# Test: replace - by dot in contrasts
gsub("-", ".", contrasts.DGE.names)

for(contrast in contrasts.DGE.names){
  contrast.dot = gsub("-", ".", contrast)
  assign(x = contrast.dot, value = topTable(ebayes.informative, coef= contrast, number=nrow(ebayes.informative), adjust.method="BH"))
  write.table(get(contrast.dot), file = file.path(topTable.dir, paste(contrast.dot, "topTable.txt", sep = ".")), quote = F, sep = "\t", row.names = T, col.names = T)
  save(list = contrast.dot, file = file.path(topTable.dir, paste(contrast.dot, "topTable.RData", sep = ".")))
}
rm(contrast.dot)



# Cleanup
rm(topTable.dir)
