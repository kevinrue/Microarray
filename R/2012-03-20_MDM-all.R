############
# Versions #
############

# R version 2.15.0 (2012-03-30)

# Warning messages:
# 1: package 'affy' was built under R version 2.15.1 
# 2: package 'BiocGenerics' was built under R version 2.15.1 
# 3: package 'Biobase' was built under R version 2.15.1
# 1: package 'bovinecdf' was built under R version 2.15.1 
# 2: package 'AnnotationDbi' was built under R version 2.15.1
# package 'arrayQualityMetrics' was built under R version 2.15.1
# 1: package 'affyPLM' was built under R version 2.15.1 
# 2: package 'gcrma' was built under R version 2.15.1 
# 3: package 'preprocessCore' was built under R version 2.15.1
# 1: package 'simpleaffy' was built under R version 2.15.1 
# 2: package 'genefilter' was built under R version 2.15.1
# 1: package 'farms' was built under R version 2.15.1 
# 2: package 'MASS' was built under R version 2.15.2
# package 'limma' was built under R version 2.15.1
# 1: package 'bovine.db' was built under R version 2.15.1 
# 2: package 'org.Hs.eg.db' was built under R version 2.15.1 
# 3: package 'RSQLite' was built under R version 2.15.2 
# package 'annotate' was built under R version 2.15.1 

##############
# Data files #
##############

#0 Used the same renamed CEL files as for the first analysis
#0 Renamed files are structured "Animal_Treatment_Timepoint.CEL"
#-1 Do not create shortcuts to the original CEL files. R wants the actual real CEL file. It does not follow shortcuts.

#############
# Libraries #
#############

library(Biobase)
library(affy) # loads the package
library(arrayQualityMetrics) # for the microarray QC
library(affyPLM) # alternative microarray QC
library(simpleaffy) #
library(farms) #
library(limma) #
library(bovine.db) #
library(annotate) #
library(puma) #
library(made4) #
library(ggplot2) #
library(PerformanceAnalytics) # for the multiple pairwise correlation neat graph

###############
# Preparation #
###############

setwd("C:\\Users\\krue\\Documents\\Kevin-Logs\\CELfiles\\Project\\ProjectData\\CEL Files detailed") # directory with CEL files
df <- data.frame(dir("."), stringsAsFactors=F) # lists all the files in the directory
names(df)[1] <- "FileName" # Renames the column
df$ID <- unlist(strsplit(df$FileName, split = ".CEL", fixed = TRUE)) # Uses the filename , withut ".CEL", as the ID
# The following block creates a column for the AnimalID, the treatment, and the timpoint
temp.list = unlist(strsplit(df$ID, split = "_", fixed = FALSE))
df$Animal <- temp.list[seq(1, length(temp.list), 3)]
df$Treatment <- temp.list[seq(2, length(temp.list), 3)]
df$TimePoint <- temp.list[seq(3, length(temp.list), 3)]

# Replace the 25HR by 24HR in the timepoint column
df[df$TimePoint=="25HR", "TimePoint"] = "24HR"

# the following block extracts the group name (treatment_timepoint) from the ID
temp.list = unlist(strsplit(df$ID, "^[0123456789]{3}R?_", fixed = FALSE, perl=TRUE))
temp.filter = temp.list != ""
df$Group <- temp.list[temp.filter]

df <- df[order(df$Treatment, df$TimePoint, df$Animal),]

setwd("C:\\Users\\krue\\Documents\\Kevin-Logs\\20130320_MDM\\R") # Directory for the current project
write.table(df, file="phenodata.txt", append=FALSE, quote=FALSE, sep="\t", eol="\n", row.names=FALSE, col.names=TRUE)

rm(temp.list, temp.filter, df) # clears the memory

###############
# Data import #
###############

setwd("C:\\Users\\krue\\Documents\\Kevin-Logs\\20130320_MDM\\R") # Directory for the current project

# read annotated file into targets information
targets = read.AnnotatedDataFrame ("phenodata.txt", header = TRUE, row.names = 2, as.is = TRUE)
#targets

# get sample names from targets
SampleNames <- sampleNames(targets)
#SampleNames

setwd("C:\\Users\\krue\\Documents\\Kevin-Logs\\CELfiles\\Project\\ProjectData\\CEL Files detailed")

# Read all CEL files into an affybatch object called rawdata
rawdata <- ReadAffy(filenames=pData(targets)$FileName, celfile.path=".", phenoData=targets, sampleNames=SampleNames)
#rawdata
#pData(rawdata)

###################
# Quality control #
###################
setwd("C:\\Users\\krue\\Documents\\Kevin-Logs\\20130320_MDM\\R") # Directory for the current project
arrayQualityMetrics(expressionset = rawdata, outdir = "MDM_raw_log", do.logtransform = TRUE)
# 700_BCG_2HR was found an outlier in 5 out of 6 tests
# 713_MAP_6HR has much higher RNA degradation and failed 2 out of 6 tests

#####################
# Samples discarded #
#####################
#1 Removed 700_BCG_2HR and 713_MAP_6HR from the phenodata file
#1 Removed all the samples for animal 700 and 713 (loss of replicates but clearer to explain to reviewers)

#################
# Data import 2 #
#################

setwd("C:\\Users\\krue\\Documents\\Kevin-Logs\\20130320_MDM\\R2") # Directory for the current project

# read annotated file into targets information
targets = read.AnnotatedDataFrame ("phenodata.txt", header = TRUE, row.names = 2, as.is = TRUE)
#targets

# get sample names from targets
SampleNames <- sampleNames(targets)
#SampleNames

setwd("C:\\Users\\krue\\Documents\\Kevin-Logs\\CELfiles\\Project\\ProjectData\\CEL Files detailed")

# Read all CEL files into an affybatch object called rawdata
rawdata <- ReadAffy(filenames=pData(targets)$FileName, celfile.path=".", phenoData=targets, sampleNames=SampleNames)
#rawdata
#pData(rawdata)

###################
# Quality control #
###################
setwd("C:\\Users\\krue\\Documents\\Kevin-Logs\\20130320_MDM\\R3") # Directory for the current project
arrayQualityMetrics(expressionset = rawdata, outdir = "MDM_raw_log", do.logtransform = TRUE)
# 721_BOVIS_24HR is barely over the threshold to be called outlier, based on Relative Log Expression. That will be fine after normalisation.

#####################
# Samples reordered #
#####################

#2 In phenodata reordered the samples by timepoint and treatments groups
#2 e.g. all controls at 0HR, all controls at 2HR, etc.
#2 Easier reading of the heatmap for distance between arrays

#################
# Data import 3 #
#################

setwd("C:\\Users\\krue\\Documents\\Kevin-Logs\\20130320_MDM\\R3") # Directory for the current project

# read annotated file into targets information
targets = read.AnnotatedDataFrame ("phenodata.txt", header = TRUE, row.names = 2, as.is = TRUE)
#targets

# get sample names from targets
SampleNames <- sampleNames(targets)
#SampleNames

setwd("C:\\Users\\krue\\Documents\\Kevin-Logs\\CELfiles\\Project\\ProjectData\\CEL Files detailed")

# Read all CEL files into an affybatch object called rawdata
rawdata <- ReadAffy(filenames=pData(targets)$FileName, celfile.path=".", phenoData=targets, sampleNames=SampleNames)
#rawdata
#pData(rawdata)

###################
# Quality control #
###################
setwd("C:\\Users\\krue\\Documents\\Kevin-Logs\\20130320_MDM\\R3") # Directory for the current project
arrayQualityMetrics(expressionset = rawdata, outdir = "MDM_raw_log", do.logtransform = TRUE)
# 721_BOVIS_24HR is barely over the threshold to be called outlier, based on Relative Log Expression. That will be fine after normalisation.


#############################
# Create boxplot of rawdata #
#############################

# NOTE: 20 samples for bovis, 20 for map
boxplot(rawdata, main="Boxplot of raw data: MAP + M.bovis", col=c(rep("pink",20),rep("lightblue",20)), las=2, cex.axis=0.53) # Create boxplot of rawdata

hist(rawdata, main="Histogram of raw data: MAP + M.bovis", col = length(rawdata):1, lty=length(rawdata):1) # Generate histogram of rawdata
legend("topright", legend=row.names(pData(rawdata)), col=length(rawdata):1, lty=length(rawdata):1, cex=0.46)

#################
# Normalisation #
#################
setwd("C:\\Users\\krue\\Documents\\Kevin-Logs\\20130320_MDM\\R3") # Directory for the current project

# Normalize the object rawdata using the FARMS normalization method. Call the normalized object farms 
farms <- qFarms(rawdata)
rm(rawdata) # don't really need it anymore

boxplot(farms, col = c(rep("pink",21),rep("lightblue",21)), las=2, cex.axis=0.53, main="Boxplot of normalised data: MAP + M.bovis") #Generate boxplot of normalised data

hist(farms, main="Histogram of normalised data: MAP + M.bovis", col = length(rawdata):1, lty=length(rawdata):1)
legend("topright", legend=row.names(pData(farms)), col=length(rawdata):1, lty=length(rawdata):1, cex=0.56)

# exports the normalised expression data
write.table(assayData(farms)$exprs, sep="\t", file="assayData-farms_exprs.txt")

#####################
# Informative probe #
#####################
setwd("C:\\Users\\krue\\Documents\\Kevin-Logs\\20130320_MDM\\R3") # Directory for the current project

# Informative probes
INIs <- INIcalls(farms) # Find informative genes from the farms normalised genes (probes?) using INIcalls

farms_informative <- getI_Eset(INIs) # List of informative probesets, as another exprSet
#farms_informative

# Save in a R environment the informative probesets and normalised data
save(farms_informative, file = "farms_informative.RData")
# Save in a text file the informative probeset and normalised data
write.table(assayData(farms_informative)$exprs, file = "farms_informative.txt")

rm(INIs) # don't really need it anymore

#########################
# Clustering of samples #
#########################
setwd("C:\\Users\\krue\\Documents\\Kevin-Logs\\20130320_MDM\\R3") # Directory for the current project

custom_col = c(rep("black",5),rep("red",5),rep("green",5),rep("blue",5),
               rep("red",5),rep("green",5),rep("blue",5),
               rep("red",5),rep("green",5),rep("blue",5),
               rep("red",5),rep("green",5),rep("blue",5))
custom_pch =c(rep(0,16), rep(1,15), rep(2,15), rep(3,15))
legend_col = c("black","red", "green","blue","red", "green","blue","red", "green","blue","red", "green","blue")
legend_pch = c(0,0,0,0,1,1,1,2,2,2,3,3,3)
## All probes (before INI calls)
pca <- prcomp(data.frame(t(assayData(farms)$exprs)), scale=T) # clustering based on expression values
#summary(pca)
#summary(pca)$importance[, 1:6]
png(filename="PCA_normalised.png", width = 600, height = 600)
plot(pca$x, col=custom_col, pch=custom_pch,
     main="PCA - Normalised Probe Sets") 
legend(x="bottomright",col=legend_col, pch=legend_pch, inset = 0.01, cex=0.60,
       legend=c("CN_0HR","CN_2HR","CN_6HR","CN_24HR","BCG_2HR","BCG_6HR","BCG_24HR",
                "MAP_2HR", "MAP_6HR", "MAP_24HR", "BOV_2HR", "BOV_6HR","BOV_24HR"))
dev.off()
# not conclusive: it takes 33 primary components to reach 80% of variance explained
# although PC1 separates very well the time points

## Informative probes only
pca <- prcomp(data.frame(t(assayData(farms_informative)$exprs)), scale=T) # clustering based on expression values
#summary(pca)
#summary(pca)$importance[, 1:6]
png(filename="PCA_informative.png", width = 600, height = 600)
plot(pca$x, col=custom_col, pch=custom_pch,
     main="PCA - Informative Probe Sets") 
legend(x="bottomright",col=legend_col, pch=legend_pch, inset = 0.01, cex=0.80,
       legend=c("CN_0HR","CN_2HR","CN_6HR","CN_24HR","BCG_2HR","BCG_6HR","BCG_24HR",
                "MAP_2HR", "MAP_6HR", "MAP_24HR", "BOV_2HR", "BOV_6HR","BOV_24HR"))
dev.off()
# slightly better: 9 primary components explain 80% of variance between arrays
# PC1 still separates very well the time points
rm(pca, custom_col, custom_pch, legend_col, legend_pch) # clears memory

###########################
# Differential expression #
###########################
setwd("C:\\Users\\krue\\Documents\\Kevin-Logs\\20130320_MDM\\R3") # Directory for the current project

groups <- unique(targets$Group) # Group animals by (treatment,timepoint)
#groups

f <- factor(targets$Group, levels = groups) # let f be the tuple (treatment,timepoint) of each array

# Design the matrix for the 65 arrays, with the columns :
# "BOVIS_2HR","BOVIS_6HR","BOVIS_24HR","MAP_2HR","MAP_6HR","MAP_24HR"
design <- model.matrix(~0 + f)
colnames(design) <- groups
#design
rm(f)

# Design the matrix listing the contrasts to compute
cntrts <- c("BOVIS_2HR-MAP_2HR","BOVIS_6HR-MAP_6HR","BOVIS_24HR-MAP_24HR",
            "BOVIS_2HR-BCG_2HR","BOVIS_6HR-BCG_6HR","BOVIS_24HR-BCG_24HR",
            "BOVIS_2HR-CN_2HR","BOVIS_6HR-CN_6HR","BOVIS_24HR-CN_24HR",
            "MAP_2HR-BCG_2HR","MAP_6HR-BCG_6HR","MAP_24HR-BCG_24HR",
            "MAP_2HR-CN_2HR","MAP_6HR-CN_6HR","MAP_24HR-CN_24HR",
            "BCG_2HR-CN_2HR","BCG_6HR-CN_6HR","BCG_24HR-CN_24HR") # Prepares the contrasts names
contrasts <- makeContrasts(contrasts=cntrts, levels=design) # Make contrast matrix showing Bovis vs Map at each time point
#contrasts

# NOTE: peculiarity for paired comparison
#unique(targets$Animal)
animal <- targets$Animal # Show Animal column in targets file and call it animal
#animal

# The following block is required to obtain the correlation value between samples from the same animal
# this value will be required to block the analysis considering non-independence between samples from each animal

##Find out what is the correlation between all arrays from the same animal. Call it array_correlation.###
# The correlation is calculated from the unfiltered list of probe sets
array_correlation <- duplicateCorrelation(farms, design=design, ndups=1, block=animal)
##Show array_correlation###
#array_correlation
#array_correlation$consensus

###Fit a linear model to the FARMS normalized data, call it farms_fit###
farms_fit = lmFit(farms, design = design, ndups=1, cor=array_correlation$consensus, block=animal)

###Apply contrasts to the linear model farms_fit###
farms_fit_contrast <- contrasts.fit(farms_fit, contrasts)

filter <- rownames(farms_fit_contrast) %in% rownames(exprs(farms_informative)) # Filter farms_fit_contrast using the informative genes found earlier (farms_informative)
filtered <- farms_fit_contrast[filter,] # long array of TRUE and FALSE whether to keep the row or not
# the above contains only the contrasted rows for the informative probesets 
# Below: Given a series of related parameter estimates and standard errors, compute moderated t-statistics, moderated F-statistic, and log-odds of differential expression by empirical Bayes shrinkage of the standard errors towards a common value.
eb_filtered <- eBayes(filtered)
#eb_filtered # Result of the differential expression statistical test, containing log odds of differential expression, t statistics, p.value, ...

# Write DEgenes to an outfile
write.table(eb_filtered, file = "DEgenes.txt", quote = FALSE, sep = "\t", eol = "\n", na = "NA", row.names = TRUE, col.names = TRUE) 
# differential expression information of informative probes stored in a file
iProbes = nrow(eb_filtered)
#number of informative probes


rm(cntrts, contrasts, animal, array_correlation, farms_fit, farms_fit_contrast, filter, filtered)

###############
# Time points #
###############
# +
######################
# bovine annotations #
######################

AnnotateTopTable <- function(df) # appends EnsEMBL and ENTREZ IDs to the tables of informative probe sets
{
  gene.symbols <- getSYMBOL(df$ID, "bovine.db")
  results <- cbind(df, gene.symbols)
  results$rankFC <- 1:length(results[,1])
  gene.entrezID <- toTable(bovineENTREZID[df$ID])
  gene.ensemblID <- toTable(bovineENSEMBL[df$ID])
  results <- merge(results, gene.entrezID, by.x="ID", by.y="probe_id", all.x=TRUE)
  results <- merge(results, gene.ensemblID, by.x="ID", by.y="probe_id", all.x=TRUE)
  # restores the ranking saved a few lines above
  results <- results[order( results$rankFC),]
  results$rankFC <- NULL
  results
}

EntrezEnsemblStat <- function(df) # counts number of NA and IDs found in EnsEMBL and ENTREZ databases (updated at each release of the bovine annotation package)
{
  Entrez.count = sum(!is.na(df$gene_id))
  Ensembl.count = sum(!is.na(df$ensembl_id))
  data.frame(cbind(Count=c(Entrez.count, Ensembl.count),
                   Percentage=c(Entrez.count/length(df$ID), Ensembl.count/length(df$ID)),
                   Total=rep(length(df$ID),2)),
             row.names=c("Entrez","Ensembl"))
}

UpsAndDowns <- function(df)# The following returns the number of up and down regulated probes in a dataset
{
  u = sum(df$logFC > 0)
  d = sum(df$logFC < 0)
  df = data.frame(rbind(Ups=u, Downs=d))
  names(df) = "bovis-MAP"
  df
}

DE <- topTable(eb_filtered, coef= "BOVIS_2HR-CN_2HR", number=iProbes, adjust.method="BH", sort.by="logFC")
DE_annotated = AnnotateTopTable(DE)
write.table(DE_annotated, file = "2Hrs_BOV_CN.txt", quote = FALSE, sep = "\t", eol = "\n", na = "NA", row.names = TRUE, col.names = TRUE)

DE <- topTable(eb_filtered, coef= "BOVIS_6HR-CN_6HR", number=iProbes, adjust.method="BH", sort.by="logFC")
DE_annotated = AnnotateTopTable(DE)
write.table(DE_annotated, file = "6Hrs_BOV_CN.txt", quote = FALSE, sep = "\t", eol = "\n", na = "NA", row.names = TRUE, col.names = TRUE)

DE <- topTable(eb_filtered, coef= "BOVIS_24HR-CN_24HR", number=iProbes, adjust.method="BH", sort.by="logFC")
DE_annotated = AnnotateTopTable(DE)
write.table(DE_annotated, file = "24Hrs_BOV_CN.txt", quote = FALSE, sep = "\t", eol = "\n", na = "NA", row.names = TRUE, col.names = TRUE)

DE <- topTable(eb_filtered, coef= "BOVIS_2HR-BCG_2HR", number=iProbes, adjust.method="BH", sort.by="logFC")
DE_annotated = AnnotateTopTable(DE)
write.table(DE_annotated, file = "2Hrs_BOV_BCG.txt", quote = FALSE, sep = "\t", eol = "\n", na = "NA", row.names = TRUE, col.names = TRUE)

DE <- topTable(eb_filtered, coef= "BOVIS_6HR-BCG_6HR", number=iProbes, adjust.method="BH", sort.by="logFC")
DE_annotated = AnnotateTopTable(DE)
write.table(DE_annotated, file = "6Hrs_BOV_BCG.txt", quote = FALSE, sep = "\t", eol = "\n", na = "NA", row.names = TRUE, col.names = TRUE)

DE <- topTable(eb_filtered, coef= "BOVIS_24HR-BCG_24HR", number=iProbes, adjust.method="BH", sort.by="logFC")
DE_annotated = AnnotateTopTable(DE)
write.table(DE_annotated, file = "24Hrs_BOV_BCG.txt", quote = FALSE, sep = "\t", eol = "\n", na = "NA", row.names = TRUE, col.names = TRUE)

DE <- topTable(eb_filtered, coef= "BOVIS_2HR-MAP_2HR", number=iProbes, adjust.method="BH", sort.by="logFC")
DE_annotated = AnnotateTopTable(DE)
write.table(DE_annotated, file = "2Hrs_BOV_MAP.txt", quote = FALSE, sep = "\t", eol = "\n", na = "NA", row.names = TRUE, col.names = TRUE)

DE <- topTable(eb_filtered, coef= "BOVIS_6HR-MAP_6HR", number=iProbes, adjust.method="BH", sort.by="logFC")
DE_annotated = AnnotateTopTable(DE)
write.table(DE_annotated, file = "6Hrs_BOV_MAP.txt", quote = FALSE, sep = "\t", eol = "\n", na = "NA", row.names = TRUE, col.names = TRUE)

DE <- topTable(eb_filtered, coef= "BOVIS_24HR-MAP_24HR", number=iProbes, adjust.method="BH", sort.by="logFC")
DE_annotated = AnnotateTopTable(DE)
write.table(DE_annotated, file = "24Hrs_BOV_MAP.txt", quote = FALSE, sep = "\t", eol = "\n", na = "NA", row.names = TRUE, col.names = TRUE)

DE <- topTable(eb_filtered, coef= "MAP_2HR-BCG_2HR", number=iProbes, adjust.method="BH", sort.by="logFC")
DE_annotated = AnnotateTopTable(DE)
write.table(DE_annotated, file = "2Hrs_MAP_BCG.txt", quote = FALSE, sep = "\t", eol = "\n", na = "NA", row.names = TRUE, col.names = TRUE)

DE <- topTable(eb_filtered, coef= "MAP_6HR-BCG_6HR", number=iProbes, adjust.method="BH", sort.by="logFC")
DE_annotated = AnnotateTopTable(DE)
write.table(DE_annotated, file = "6Hrs_MAP_BCG.txt", quote = FALSE, sep = "\t", eol = "\n", na = "NA", row.names = TRUE, col.names = TRUE)

DE <- topTable(eb_filtered, coef= "MAP_24HR-BCG_24HR", number=iProbes, adjust.method="BH", sort.by="logFC")
DE_annotated = AnnotateTopTable(DE)
write.table(DE_annotated, file = "24Hrs_MAP_BCG.txt", quote = FALSE, sep = "\t", eol = "\n", na = "NA", row.names = TRUE, col.names = TRUE)

DE <- topTable(eb_filtered, coef= "MAP_2HR-CN_2HR", number=iProbes, adjust.method="BH", sort.by="logFC")
DE_annotated = AnnotateTopTable(DE)
write.table(DE_annotated, file = "2Hrs_MAP_CN.txt", quote = FALSE, sep = "\t", eol = "\n", na = "NA", row.names = TRUE, col.names = TRUE)

DE <- topTable(eb_filtered, coef= "MAP_6HR-CN_6HR", number=iProbes, adjust.method="BH", sort.by="logFC")
DE_annotated = AnnotateTopTable(DE)
write.table(DE_annotated, file = "6Hrs_MAP_CN.txt", quote = FALSE, sep = "\t", eol = "\n", na = "NA", row.names = TRUE, col.names = TRUE)

DE <- topTable(eb_filtered, coef= "MAP_24HR-CN_24HR", number=iProbes, adjust.method="BH", sort.by="logFC")
DE_annotated = AnnotateTopTable(DE)
write.table(DE_annotated, file = "24Hrs_MAP_CN.txt", quote = FALSE, sep = "\t", eol = "\n", na = "NA", row.names = TRUE, col.names = TRUE)

DE <- topTable(eb_filtered, coef= "BCG_2HR-CN_2HR", number=iProbes, adjust.method="BH", sort.by="logFC")
DE_annotated = AnnotateTopTable(DE)
write.table(DE_annotated, file = "2Hrs_BCG_CN.txt", quote = FALSE, sep = "\t", eol = "\n", na = "NA", row.names = TRUE, col.names = TRUE)

DE <- topTable(eb_filtered, coef= "BCG_6HR-CN_6HR", number=iProbes, adjust.method="BH", sort.by="logFC")
DE_annotated = AnnotateTopTable(DE)
write.table(DE_annotated, file = "6Hrs_BCG_CN.txt", quote = FALSE, sep = "\t", eol = "\n", na = "NA", row.names = TRUE, col.names = TRUE)

DE <- topTable(eb_filtered, coef= "BCG_24HR-CN_24HR", number=iProbes, adjust.method="BH", sort.by="logFC")
DE_annotated = AnnotateTopTable(DE)
write.table(DE_annotated, file = "24Hrs_BCG_CN.txt", quote = FALSE, sep = "\t", eol = "\n", na = "NA", row.names = TRUE, col.names = TRUE)

rm(DE)

####################
# Packages updated #
####################

#3 The annotation package and Bioconductor in general was not up to date
#3 Updated and repeated the same analysis as above in R4

# R version 2.15.3 (2013-03-01) -- "Security Blanket"
# Only warning message:
# "No methods found in "Biobase" for requests: geneNames"
# Not sure whether this is important




#####################################################################################################
#####################################################################################################



################### This is the last, definitive way of importing data that I've done.
################### It only includes 5 animals passing QC for all samples.
## Data import 4 ##
###################
###################
setwd("C:\\Users\\krue\\Documents\\Kevin-Logs\\20130320_MDM\\R4_QC-filtered") # Directory for the current project

# read annotated file into targets information
targets = read.AnnotatedDataFrame ("phenodata.txt", header = TRUE, row.names = 2, as.is = TRUE)
#targets

# get sample names from targets
SampleNames <- sampleNames(targets)
#SampleNames

setwd("C:\\Users\\krue\\Documents\\Kevin-Logs\\CELfiles\\Project\\ProjectData\\CEL Files detailed")

# Read all CEL files into an affybatch object called rawdata
rawdata <- ReadAffy(filenames=pData(targets)$FileName, celfile.path=".", phenoData=targets, sampleNames=SampleNames)
#rawdata
#pData(rawdata)

###################
# Quality control #
###################
setwd("C:\\Users\\krue\\Documents\\Kevin-Logs\\20130320_MDM\\R4_QC-filtered") # Directory for the current project
arrayQualityMetrics(expressionset = rawdata, outdir = "MDM_raw_log", do.logtransform = TRUE)
# 721_BOVIS_24HR is barely over the threshold to be called outlier, based on Relative Log Expression. That will be fine after normalisation.

#############################
# Create boxplot of rawdata #
#############################
# Repetitive of the arrayQualityMetrics package

# NOTE: 20 samples for bovis, 20 for map
boxplot(rawdata, main="Boxplot of raw data: MDM experiment", col=c(rep("pink",20),rep("lightblue",20)), las=2, cex.axis=0.53) # Create boxplot of rawdata

hist(rawdata, main="Histogram of raw data: MDM experiment", col = length(rawdata):1, lty=length(rawdata):1) # Generate histogram of rawdata
legend("topright", legend=row.names(pData(rawdata)), col=length(rawdata):1, lty=length(rawdata):1, cex=0.46)

#################
# Normalisation #
#################
setwd("C:\\Users\\krue\\Documents\\Kevin-Logs\\20130320_MDM\\R4_QC-filtered") # Directory for the current project

# Normalize the object rawdata using the FARMS normalization method. Call the normalized object farms 
farms <- qFarms(rawdata)
rm(rawdata) # don't really need it anymore

#Generate boxplot of normalised data
boxplot(farms, names=SampleNames, col = c(rep("pink",65)),las=2, cex.axis=0.53, main="Boxplot of normalised data: MAP + M.bovis")

hist(farms, main="Histogram of normalised data: MDM experiment", col = length(SampleNames):1, lty=length(SampleNames):1)
legend("topright", legend=row.names(pData(farms)), col=length(SampleNames):1, lty=length(SampleNames):1, cex=0.60)

# exports the normalised expression data
write.table(assayData(farms)$exprs, sep="\t", file="assayData-farms_exprs.txt")

#####################
# Informative probe #
#####################
setwd("C:\\Users\\krue\\Documents\\Kevin-Logs\\20130320_MDM\\R4") # Directory for the current project

# Informative probes
INIs <- INIcalls(farms) # Find informative genes from the farms normalised genes (probes?) using INIcalls

farms_informative <- getI_Eset(INIs) # List of informative probesets, as another exprSet
#farms_informative

# Save in a R environment the informative probesets and normalised data
save(farms_informative, file = "farms_informative.RData")
# Save in a text file the informative probeset and normalised data
write.table(assayData(farms_informative)$exprs, file = "farms_informative.txt")

rm(INIs) # don't really need it anymore

#########################
# Clustering of samples #
#########################
setwd("C:\\Users\\krue\\Documents\\Kevin-Logs\\20130320_MDM\\R4_QC-filtered") # Directory for the current project

custom_col = c(rep("black",5),rep("red",5),rep("green",5),rep("blue",5),
               rep("red",5),rep("green",5),rep("blue",5),
               rep("red",5),rep("green",5),rep("blue",5),
               rep("red",5),rep("green",5),rep("blue",5))
custom_pch =c(rep(0,16), rep(1,15), rep(2,15), rep(3,15))
#custom_pch =c(rep(4,16), rep(5,15), rep(6,15), rep(7,15))
#custom_pch =c(rep(8,16), rep(9,15), rep(10,15), rep(11,15))
#custom_pch =c(rep(12,16), rep(13,15), rep(14,15), rep(15,15))

#custom_cex =c(rep(0.5,16), rep(1,15), rep(1.5,15), rep(2,15))
legend_col = c("black","red", "green","blue","red", "green","blue","red", "green","blue","red", "green","blue")
legend_pch = c(0,0,0,0,1,1,1,2,2,2,3,3,3)

##
## All probes (before INI calls) ##
#pca <- prcomp(data.frame(t(assayData(farms)$exprs)), scale=T) # clustering based on expression values
#summary(pca)
#summary(pca)$importance[, 1:6]
#png(filename="PCA_normalised.png", width = 600, height = 600)
#plot(pca$x, col=custom_col, pch=custom_pch,
#     main="PCA - Normalised Probe Sets") 
#legend(x="bottomright",col=legend_col, pch=legend_pch, inset = 0.01, cex=0.60,
#       legend=c("CN_0HR","CN_2HR","CN_6HR","CN_24HR","BCG_2HR","BCG_6HR","BCG_24HR",
                "MAP_2HR", "MAP_6HR", "MAP_24HR", "BOV_2HR", "BOV_6HR","BOV_24HR"))
#dev.off()
# not conclusive: it takes 33 primary components to reach 80% of variance explained
# although PC1 separates very well the time points

##
## Informative probes only ##
pca <- prcomp(data.frame(t(assayData(farms_informative)$exprs)), scale=T) # clustering based on expression values

#summary(pca)
#summary(pca)$importance[, 1:6]
png(filename="PCA_informative.png", width = 600, height = 600)
plot(pca$x, col=custom_col, pch=custom_pch,
     main="PCA - Informative Probe Sets") 
legend(x="bottomright",col=legend_col, pch=legend_pch, inset = 0.01, cex=0.80,
       legend=c("CN_0HR","CN_2HR","CN_6HR","CN_24HR","BCG_2HR","BCG_6HR","BCG_24HR",
                "MAP_2HR", "MAP_6HR", "MAP_24HR", "BOV_2HR", "BOV_6HR","BOV_24HR"))
dev.off()
# slightly better: 9 primary components explain 80% of variance between arrays
# PC1 still separates very well the time points

rm(pca, custom_col, custom_pch, legend_col, legend_pch) # clears memory


#####
# 
####

# color coding
pca_coord = data.frame(pca$x)
# time coloring
pca_coord$time.col = NA
pca_coord$time.col[grep("0HR", rownames(pca_coord))] = "black"
pca_coord$time.col[grep("2HR", rownames(pca_coord))] = "red"
pca_coord$time.col[grep("6HR", rownames(pca_coord))] = "green"
pca_coord$time.col[grep("24HR", rownames(pca_coord))] = "blue"
pca_coord$time.col[grep("25HR", rownames(pca_coord))] = "blue" # names 25HR are considered 24HR
# treatment coloring
pca_coord$treat.col = NA
pca_coord$treat.col[grep("CN", rownames(pca_coord))] = "black"
pca_coord$treat.col[grep("BCG", rownames(pca_coord))] = "red"
pca_coord$treat.col[grep("MAP", rownames(pca_coord))] = "green"
pca_coord$treat.col[grep("BOVIS", rownames(pca_coord))] = "blue"
# time shaping
pca_coord$time.pch = NA
pca_coord$time.pch[grep("0HR", rownames(pca_coord))] = 1
pca_coord$time.pch[grep("2HR", rownames(pca_coord))] = 3
pca_coord$time.pch[grep("6HR", rownames(pca_coord))] = 10
pca_coord$time.pch[grep("24HR", rownames(pca_coord))] = 11
pca_coord$time.pch[grep("25HR", rownames(pca_coord))] = 11 # names 25HR are considered 24HR

# Tried PCs
# 1-2, 1-3, 2-3, 3-4, 2-4
# Problem is I can only assess by eye.
# I describe a possible systematic way of doing it, in my report

plot(pca_coord$PC4, pca_coord$PC3, col=pca_coord$treat.col, pch=pca_coord$time.pch,
     main="PCA - Informative Probe Sets") 
legend(x="topright",col=c(rep("black",4), rep("red",3), rep("green",3), rep("blue",3)),
       pch=c(1,3,10,11,3,10,11,3,10,11), inset = 0.01, cex=0.80,
       legend=c("CN_0HR","CN_2HR","CN_6HR","CN_24HR","BCG_2HR","BCG_6HR","BCG_24HR",
                "MAP_2HR", "MAP_6HR", "MAP_24HR", "BOV_2HR", "BOV_6HR","BOV_24HR"))

##########
# MADE 4 #
##########
# Paul described it as a tool which can neatly represent pre-defined clusters
# in expression-based PCA plots.

# I may use it for a more proper graphical representation of my clusters.

###########################
# Differential expression #
###########################
setwd("C:\\Users\\krue\\Documents\\Kevin-Logs\\20130320_MDM\\R4") # Directory for the current project

groups <- unique(targets$Group) # Group animals by (treatment,timepoint)
#groups

f <- factor(targets$Group, levels = groups) # let f be the tuple (treatment,timepoint) of each array

# Design the matrix for the 65 arrays, with the columns :
# "BOVIS_2HR","BOVIS_6HR","BOVIS_24HR","MAP_2HR","MAP_6HR","MAP_24HR"
design <- model.matrix(~0 + f)
colnames(design) <- groups
#design
rm(f)

# Design the matrix listing the contrasts to compute
cntrts <- c("BOVIS_2HR-MAP_2HR","BOVIS_6HR-MAP_6HR","BOVIS_24HR-MAP_24HR",
            "BOVIS_2HR-BCG_2HR","BOVIS_6HR-BCG_6HR","BOVIS_24HR-BCG_24HR",
            "BOVIS_2HR-CN_2HR","BOVIS_6HR-CN_6HR","BOVIS_24HR-CN_24HR",
            "MAP_2HR-BCG_2HR","MAP_6HR-BCG_6HR","MAP_24HR-BCG_24HR",
            "MAP_2HR-CN_2HR","MAP_6HR-CN_6HR","MAP_24HR-CN_24HR",
            "BCG_2HR-CN_2HR","BCG_6HR-CN_6HR","BCG_24HR-CN_24HR") # Prepares the contrasts names
contrasts <- makeContrasts(contrasts=cntrts, levels=design) # Make contrast matrix showing Bovis vs Map at each time point
#contrasts

# NOTE: peculiarity for paired comparison
#unique(targets$Animal)
animal <- targets$Animal # Show Animal column in targets file and call it animal
#animal

# The following block is required to obtain the correlation value between samples from the same animal
# this value will be required to block the analysis considering non-independence between samples from each animal

##Find out what is the correlation between all arrays from the same animal. Call it array_correlation.###
# The correlation is calculated from the unfiltered list of probe sets
array_correlation <- duplicateCorrelation(farms, design=design, ndups=1, block=animal)
##Show array_correlation###
#array_correlation
#array_correlation$consensus

###Fit a linear model to the FARMS normalized data, call it farms_fit###
farms_fit = lmFit(farms, design = design, ndups=1, cor=array_correlation$consensus, block=animal)

###Apply contrasts to the linear model farms_fit###
farms_fit_contrast <- contrasts.fit(farms_fit, contrasts)

filter <- rownames(farms_fit_contrast) %in% rownames(exprs(farms_informative)) # Filter farms_fit_contrast using the informative genes found earlier (farms_informative)
filtered <- farms_fit_contrast[filter,] # long array of TRUE and FALSE whether to keep the row or not
# the above contains only the contrasted rows for the informative probesets 
# Below: Given a series of related parameter estimates and standard errors, compute moderated t-statistics, moderated F-statistic, and log-odds of differential expression by empirical Bayes shrinkage of the standard errors towards a common value.
eb_filtered <- eBayes(filtered)
#eb_filtered # Result of the differential expression statistical test, containing log odds of differential expression, t statistics, p.value, ...

# Write DEgenes to an outfile
write.table(eb_filtered, file = "DEgenes - eb_filtered.txt", quote = FALSE, sep = "\t", eol = "\n", na = "NA", row.names = TRUE, col.names = TRUE) 
# differential expression information of informative probes stored in a file
iProbes = nrow(eb_filtered)
#number of informative probes


rm(contrasts, animal, array_correlation, farms_fit, farms_fit_contrast, filter, filtered)

###############
# Time points #
###############
# +
######################
# bovine annotations #
######################

AnnotateTopTable <- function(df) # appends EnsEMBL and ENTREZ IDs to the tables of informative probe sets
{
  gene.symbols <- getSYMBOL(df$ID, "bovine.db")
  results <- cbind(df, gene.symbols)
  results$rank <- 1:length(results[,1])
  gene.entrezID <- toTable(bovineENTREZID[df$ID])
  gene.ensemblID <- toTable(bovineENSEMBL[df$ID])
  results <- merge(results, gene.entrezID, by.x="ID", by.y="probe_id", all.x=TRUE)
  results <- merge(results, gene.ensemblID, by.x="ID", by.y="probe_id", all.x=TRUE)
  # restores the ranking saved a few lines above
  results <- results[order( results$rank),]
  results$rank <- NULL
  results
}

EntrezEnsemblStat <- function(df) # counts number of NA and IDs found in EnsEMBL and ENTREZ databases (updated at each release of the bovine annotation package)
{
  Entrez.count = sum(!is.na(df$gene_id))
  Ensembl.count = sum(!is.na(df$ensembl_id))
  data.frame(cbind(Count=c(Entrez.count, Ensembl.count),
                   Percentage=c(Entrez.count/length(df$ID), Ensembl.count/length(df$ID)),
                   Total=rep(length(df$ID),2)),
             row.names=c("Entrez","Ensembl"))
}

UpsAndDowns <- function(df, cntr, out)# The following returns the number of up and down regulated probes in a dataset
{
  # analysis based on classic p value
  up.FC0p05 = sum(df$logFC > 0 & df$adj.P.Val < 0.05)
  down.FC0p05 = sum(df$logFC < 0 & df$adj.P.Val < 0.05)
  DE.FC0p05 = sum(df$adj.P.Val < 0.05)
  nDE.FC0p05 = sum(df$adj.P.Val > 0.05)
  # increased strigency on Log Fold Change
  up.FC5p05 = sum(df$logFC > 0.5 & df$adj.P.Val < 0.05)
  down.FC5p05 = sum(df$logFC < -0.5 & df$adj.P.Val < 0.05)
  DE.FC5p05 = sum(abs(df$logFC) > 0.5 & df$adj.P.Val < 0.05)
  nDE.FC5p05 = sum(abs(df$logFC < 0.5) | df$adj.P.Val > 0.05)
  # increased strigency on p  value
  up.FC0p01 = sum(df$logFC > 0 & df$adj.P.Val < 0.01)
  down.FC0p01 = sum(df$logFC < 0 & df$adj.P.Val < 0.01)
  DE.FC0p01 = sum(df$adj.P.Val < 0.01)
  nDE.FC0p01 = sum(df$adj.P.Val > 0.01)
  # put data in data frame 
  out = rbind(out, c(up.FC0p05,down.FC0p05,DE.FC0p05,nDE.FC0p05,
                     up.FC5p05,down.FC5p05,DE.FC5p05,nDE.FC5p05,
                     up.FC0p01,down.FC0p01,DE.FC0p01,nDE.FC0p01))
}

# For each contrast:
## get the full list of informative probes, sorted by P-value
## annotate the ordered list of probes with EnsEMBL and ENTREZ IDs
## add an entry in a summary table, of genes up and down (p-value or 0.05 and 0.01 should be enough)

# cntrts contains the list of contrast names. Use it to loop through the results

DE.df = data.frame()
for(cntr in cntrts){
  DE = topTable(eb_filtered, coef= cntr, number=iProbes, adjust.method="BH", sort.by="P")
  DE.df = UpsAndDowns(DE, cntr, DE.df)
  DE_annotated = AnnotateTopTable(DE)
  write.table(DE_annotated, file = paste(cntr,".txt", sep=""), quote = FALSE, sep = "\t", eol = "\n", na = "NA", row.names = FALSE, col.names = TRUE)
}
rm(DE, cntr, DE_annotated)
rownames(DE.df) = cntrts
colnames(DE.df) = c("up.FC0p05","down.FC0p05","DE.FC0p05","nDE.FC0p05",
                    "up.FC5p05","down.FC5p05","DE.FC5p05","nDE.FC5p05",
                    "up.FC0p01","down.FC0p01","DE.FC0p01","nDE.FC0p01")
DE.df
write.table(x=DE.df, file="DEgenes - comparison-threshold.txt", quote=FALSE, sep="\t", eol = "\n", na = "NA", row.names = TRUE, col.names = TRUE)
rm(DE.df)

################
# PUMA package #
################
## library(PUMA)
# Paul suggested that another linear model considering time and treatment as separate factors may call DE genes better
# and especially help us find differences between MAP and BCG

# prepare the data for Puma
farms2 = farms  # duplicate the data because the puma functions are more restrictive about the data layout
farms2 = farms2[,which(farms2$TimePoint != "0HR")] # the zero time point prevent comparison across treatment and time points
pData(farms2) = pData(farms)[pData(farms)$TimePoint != "0HR",3:4] # the extra columns are confusing for eBayes which does not need them
# the actual analysis
design = createDesignMatrix(farms2)
fit = lmFit(farms2, design)
contM = createContrastMatrix(farms2)
fitC = contrasts.fit(fit,contM)
fitC = eBayes(fitC)
# Get the number DE calls at q-value 0.05
esClas = decideTests(fitC,p.value=0.05,adjust.method="BH")
# Get the total count of DE genes for each contrast
#colSums(abs(esClas))
write.table(x=colSums(abs(esClas)), file="PUMA Multiple contrasts DE calls.txt", append=FALSE, quote=FALSE, sep="\t", eol="\n", dec=".", col.names="DE.counts")

rm(farms2, design, fit, contM, fitC)

#####################
# log2FC 0 / Q 0.05 #
#####################
# I decided to run a first analysis at the less stringent threshold
# Questions of interest:
## Union of probe sets DE in at least one of the comparison? 
## Gene ontology of the above list?
## Union of probe sets specific to given conditions (e.g. union of all time points for BOVIS-CN, minus probes found DE in MAP-CN and BCG-CN)

# Summary of which probe sets is DE in which condition (useful for subsetting and Venn diagram)
# Note:  it will only conmsider one threshold: Q-value 0.05 and no threshold on fold-change
# probe name
# TRUE/FALSE for each contrast
# sum of TRUE
# gene.symbol
# ENTREZ+ENSEMBL IDs
DE.count = data.frame(row.names=1:iProbes) # initialises the data frame
DE.count$ID = topTable(eb_filtered, coef=1, number=iProbes, adjust.method="BH")$ID # gets the IDs of the informative probe sets
# Progressively merges the ID table with information about each contrast
for(cntr in cntrts){
  topT = topTable(eb_filtered, coef= cntr, number=iProbes, adjust.method="BH")
  topT$test = topT$adj.P.Val < 0.05
  names(topT)[names(topT)=="test"] <- cntr
  DE.count = merge(x=DE.count, y=topT[c("ID",cntr)], by="ID", sort=TRUE)
}
# clear memory
rm(topT)
# Summary stat: number of contrasts where each probe set was DE
DE.count$ContrastSum = rowSums(DE.count[names(DE.count)[grep(pattern="-", x=names(DE.count))]])
# Gene Annotation
iDE_count = AnnotateTopTable(DE.df) # 11952 rows, with duplicates from the 11842 informative probe sets
# write it to a file
write.table(x=iDE_count, file="Informative_probes_DE_contrasts.csv", append=F, quote=F, sep=",", eol="\n", dec=".", row.names=F)


## Union of probe sets DE in at least one of the comparison? 
# Be careful for count. Work with the non-annotated data frames
# as the annotated ones have dplicated rows due to multiple IDs mapped to soem probe sets
nrow(DE.count[DE.count$ContrastSum > 0,]) # 8113 probe sets are found DE at 1+ time points / out of 11842 informative probe sets / out of 24128 probe sets on the array
Positive_DE_count = AnnotateTopTable(DE.count[DE.count$ContrastSum > 0,]) # As I wrote above, this data set has 8193 rows, now including duplicates
write.table(x=Positive_DE_count, file="Positive_probes_DE_contrasts.csv", append=F, quote=F, sep=",", eol="\n", dec=".", row.names=F)


## Specific conditions: MAP/CN vs BCP/CN
count_subsets = function(df, cntr1, cntr2)
{
  # given two conditions, returns how many samples are DE in each 4 possible cases (none, 1, 2, 1 and 2)
  result = matrix(nrow=2,ncol=2,dimnames=list(c(cntr1, paste("-",cntr1, sep="")),c(cntr2, paste("-",cntr2, sep=""))))
  result[1,1] = nrow(df[df[cntr1] & df[cntr2],]) # True/True
  result[1,2] = nrow(df[df[cntr1] & !df[cntr2],]) # True/False
  result[2,1] = nrow(df[!df[cntr1] & df[cntr2],]) # False/True
  result[2,2] = nrow(df[!df[cntr1] & !df[cntr2],]) # False/False
  # margins
  addmargins(result)
}
# We want to know how many of the DE genes overlap (or not) in two given contrasts
# table counts
count_subsets(DE.count, "MAP_2HR-CN_2HR", "BCG_2HR-CN_2HR")
count_subsets(DE.count, "MAP_6HR-CN_6HR", "BCG_6HR-CN_6HR")
count_subsets(DE.count, "MAP_24HR-CN_24HR", "BCG_24HR-CN_24HR")
# the above is deprecated, see below
# which is a graphical representation of it without need of the homemade function (does require clustering of the probe sets)
vennDiagram(DE.count[,c("MAP_2HR-CN_2HR", "BCG_2HR-CN_2HR")])
vennDiagram(DE.count[,c("MAP_6HR-CN_6HR", "BCG_6HR-CN_6HR")])
vennDiagram(DE.count[,c("MAP_24HR-CN_24HR", "BCG_24HR-CN_24HR")])

## Genes unique to BCG/CN at 6 hours (not found in MAP/CN)
# At 6 hpi, there is a particularly high number of DE genes in BCG/CN
# compared to MAP/CN, while none of these show up in MAP/BCG contrast.
# What is the expression profile of these genes, which makes them being 
# called DE in BCG/CN but neither MAP/CN nor BCG/MAP
tmp = DE.count[!DE.count["MAP_6HR-CN_6HR"] & DE.count["BCG_6HR-CN_6HR"],]
tmp_ID = tmp$ID

plot_expression = function(probe, gene)
{
  # Plots the normalised expression values for each condition for a given
  # probe
  # creates a temporary copy of the expression data, annotated by condition
  tmp_df = data.frame(t(assayData(farms_informative)$exprs)) #
  tmp_df$group = NA #
  tmp_df$group[grep("CN_0HR",rownames(tmp_df))] = "CN_0HR" #
  tmp_df$group[grep("CN_2HR",rownames(tmp_df))] = "CN_2HR" #
  tmp_df$group[grep("CN_6HR",rownames(tmp_df))] = "CN_6HR" #
  tmp_df$group[grep("CN_24HR",rownames(tmp_df))] = "CN_24HR" #
  tmp_df$group[grep("CN_25HR",rownames(tmp_df))] = "CN_24HR" #
  
  tmp_df$group[grep("BCG_2HR",rownames(tmp_df))] = "BCG_2HR" #
  tmp_df$group[grep("BCG_6HR",rownames(tmp_df))] = "BCG_6HR" #
  tmp_df$group[grep("BCG_24HR",rownames(tmp_df))] = "BCG_24HR" #
  tmp_df$group[grep("BCG_25HR",rownames(tmp_df))] = "BCG_24HR" #
  
  tmp_df$group[grep("MAP_2HR",rownames(tmp_df))] = "MAP_2HR" #
  tmp_df$group[grep("MAP_6HR",rownames(tmp_df))] = "MAP_6HR" #
  tmp_df$group[grep("MAP_24HR",rownames(tmp_df))] = "MAP_24HR" #
  tmp_df$group[grep("MAP_25HR",rownames(tmp_df))] = "MAP_24HR" #
  
  tmp_df$group[grep("BOVIS_2HR",rownames(tmp_df))] = "BOVIS_2HR" #
  tmp_df$group[grep("BOVIS_6HR",rownames(tmp_df))] = "BOVIS_6HR" #
  tmp_df$group[grep("BOVIS_24HR",rownames(tmp_df))] = "BOVIS_24HR" #
  tmp_df$group[grep("BOVIS_25HR",rownames(tmp_df))] = "BOVIS_24HR" #
  
  #tmp_df$group = as.factor(tmp_df$group) #
  tmp_df$group = factor(tmp_df$group, groups) # alternative way, but "groups" has to be accurate!
  boxplot(tmp_df[,probe]~group, data=tmp_df, cex.axis=0.78, las=3, main=paste(gene, "=", probe),
          col=c(rep("green",4), rep("blue",3), rep("brown",3), rep("red",3)))
}

# Actually it is just faster by hand:
par(mfrow=c(5,5), ask=TRUE)
# groups samples by conditions
tmp_df = data.frame(t(assayData(farms_informative)$exprs)) #
tmp_df$group = NA #
tmp_df$group[grep("CN_0HR",rownames(tmp_df))] = "CN_0HR" #
tmp_df$group[grep("CN_2HR",rownames(tmp_df))] = "CN_2HR" #
tmp_df$group[grep("CN_6HR",rownames(tmp_df))] = "CN_6HR" #
tmp_df$group[grep("CN_24HR",rownames(tmp_df))] = "CN_24HR" #
tmp_df$group[grep("CN_25HR",rownames(tmp_df))] = "CN_24HR" #

tmp_df$group[grep("BCG_2HR",rownames(tmp_df))] = "BCG_2HR" #
tmp_df$group[grep("BCG_6HR",rownames(tmp_df))] = "BCG_6HR" #
tmp_df$group[grep("BCG_24HR",rownames(tmp_df))] = "BCG_24HR" #
tmp_df$group[grep("BCG_25HR",rownames(tmp_df))] = "BCG_24HR" #

tmp_df$group[grep("MAP_2HR",rownames(tmp_df))] = "MAP_2HR" #
tmp_df$group[grep("MAP_6HR",rownames(tmp_df))] = "MAP_6HR" #
tmp_df$group[grep("MAP_24HR",rownames(tmp_df))] = "MAP_24HR" #
tmp_df$group[grep("MAP_25HR",rownames(tmp_df))] = "MAP_24HR" #

tmp_df$group[grep("BOVIS_2HR",rownames(tmp_df))] = "BOVIS_2HR" #
tmp_df$group[grep("BOVIS_6HR",rownames(tmp_df))] = "BOVIS_6HR" #
tmp_df$group[grep("BOVIS_24HR",rownames(tmp_df))] = "BOVIS_24HR" #
tmp_df$group[grep("BOVIS_25HR",rownames(tmp_df))] = "BOVIS_24HR" #

tmp_df$group = as.factor(tmp_df$group) #
tmp_df$group = factor(tmp_df$group, groups) #
# box-plots the expression values by the above-defined group
for(ID in tmp_ID[1:50])
{
  boxplot(tmp_df[,ID]~group, data=tmp_df, cex.axis=0.7, main=ID,
          col=c(rep("green",4), rep("blue",3), rep("brown",3), rep("red",3)))
}
par(mfrow=c(1,1), ask=FALSE)
rm(tmp_df)

# For demonstration purpose, random selection of 9 graphs (3x3 figure)
par(mfrow=c(3,3))
test = sample(1:420, 9)
png(filename="Sample_9_Expression-BCG-not-MAP-vs-CN.png", width=1700, height=1000)
par(mfrow=c(3,3))
for(ID in tmp_ID[test])
{
  boxplot(tmp_df[,ID]~group, data=tmp_df, cex.axis=1, main=ID, las=2,
          col=c(rep("green",4), rep("blue",3), rep("brown",3), rep("red",3)))
}
dev.off()
rm(tmp_ID, tmp_df, ID)


## Venn diagram of BOV-CN, MAP-CN and BCG-CN
vennDiagram(DE.count[,c("MAP_2HR-CN_2HR", "BCG_2HR-CN_2HR", "BOVIS_2HR-CN_2HR")])
vennDiagram(DE.count[,c("MAP_6HR-CN_6HR", "BCG_6HR-CN_6HR", "BOVIS_6HR-CN_6HR")])
vennDiagram(DE.count[,c("MAP_24HR-CN_24HR", "BCG_24HR-CN_24HR", "BOVIS_24HR-CN_24HR")])


#######
# IPA #
#######
setwd("C:\\Users\\krue\\Documents\\Kevin-Logs\\20130320_MDM\\IPA_datasets-test")

# Testing purpose. 
# Does IPA obtain similar results from a dataset containing all informative probe sets
# compared to one pre-filtered for informative probes with q < 0.05 ?
# Testing on BCG-CN at 6HR (reasonable number of DE genes)

## Dataset of all informative probe sets
setwd("C:\\Users\\krue\\Documents\\Kevin-Logs\\20130320_MDM\\IPA_datasets-test")
DE = topTable(eb_filtered, coef= "BCG_6HR-CN_6HR", number=iProbes, adjust.method="BH", sort.by="P")
write.table(DE, file = paste("BCG_6HR-CN_6HR",".txt", sep=""), quote = FALSE, sep = "\t", eol = "\n", na = "NA", row.names = FALSE, col.names = TRUE)
## Dataset of pre-filtered informative probe sets
DE = topTable(eb_filtered, coef= "BCG_6HR-CN_6HR", number=iProbes, adjust.method="BH", sort.by="P", p.value=0.05)
write.table(DE, file = paste("BCG_6HR-CN_6HR_p05",".txt", sep=""), quote = FALSE, sep = "\t", eol = "\n", na = "NA", row.names = FALSE, col.names = TRUE)

# Actual project
# Generated all non-redundant dataset (no gene annotation, just probe ID)
# Therefore same number of rows as informative probe sets
# The filtering will be done ib IPA at the moment of the batch analysis of all datasets 
# adj p-value of 0.05
# MAP_BCG at 2 and 6 hours don not have any DE genes, thus no IPA analysis possible for those two data sets. --> 16 datasets processed
# The 18 datasets were still uploaded for comprehensiveness
setwd("C:\\Users\\krue\\Documents\\Kevin-Logs\\20130320_MDM\\IPA_contrasts_informative")
for(cntr in cntrts){
  DE = topTable(eb_filtered, coef= cntr, number=iProbes, adjust.method="BH", sort.by="P")
  write.table(DE, file = paste(cntr,".txt", sep=""), quote = FALSE, sep = "\t", eol = "\n", na = "NA", row.names = FALSE, col.names = TRUE)
}
rm(DE, cntr)

###########################
# Specific probe boxplots #
###########################
setwd("C:\\Users\\krue\\Documents\\Kevin-Logs\\20130320_MDM\\Expression_boxplots")
# IL1B
png(filename="Bt.4856.1.S1_at_IL1B.png", width=700, 700)
plot_expression("Bt.4856.1.S1_at", "IL1B")
dev.off()

# As a faster function
boxplot_expression = function(probe, gene)
{
  png(filename=paste(paste(probe, gene, sep="_"), "png", sep="."), width=700, 700)
  plot_expression(probe, gene)
  dev.off()
}

boxplot_expression("Bt.4856.1.S1_at", "IL1B")
boxplot_expression("Bt.4856.1.S2_at", "IL1B")
# All the other probe sets with pattern "*IFN*"
boxplot_expression("Bt.416.1.S1_at", "IFNAC")
boxplot_expression("BtAffx.1.7.S1_at", "IFNB1")
boxplot_expression("Bt.188.1.S1_at", "IFNG")
boxplot_expression("Bt.4557.1.S1_at", "IFNAR1")
boxplot_expression("Bt.5508.1.S1_at", "IFNAR2")
boxplot_expression("Bt.4251.1.S1_at", "IFNGR2")
boxplot_expression("Bt.4251.2.A1_at", "IFNGR2")


plot_expression("Bt.4856.1.S1_at", "IL1B")
plot_expression("Bt.26983.1.S1_at","CASP1")
plot_expression("Bt.191.1.S1_at","IL1A")
plot_expression("Bt.191.1.S2_at","IL1A")
plot_expression("Bt.16966.1.S1_at","CXCL10")
# Type 1 Interferon pathway
plot_expression("Bt.416.1.S1_at", "IFNAC")
plot_expression("BtAffx.1.7.S1_at", "IFNB1")
plot_expression("Bt.16966.1.S1_at","CXCL10")
# Type 2 interferon
plot_expression("Bt.188.1.S1_at", "IFNG")


# I need a function which tells whether a given probe is DE in each contrast
# Luckily, I already saved this information in a table (DE.count)
# I just need the function to retrieve the corresponding row
get_DE_counts = function(probe, gene)
{
  # gene is just here to force typing it, and thus recording which probe is which gene
  t(DE.count[DE.count$ID == probe,])
}

# IL-1 signalling
get_DE_counts("Bt.4856.1.S1_at","IL1B")
get_DE_counts("Bt.191.1.S1_at","IL1A")

# Type 1 interferon signalling
get_DE_counts("BtAffx.1.7.S1_at", "IFNB1")
get_DE_counts("Bt.416.1.S1_at", "IFNAC")
get_DE_counts("Bt.16966.1.S1_at","CXCL10")

# Type 2 interferon
get_DE_counts("Bt.188.1.S1_at", "IFNG")

# Caspase-1
get_DE_counts("Bt.26983.1.S1_at","CASP1")

# IFN-inducible genes (Novikov, 2011)
##get_DE_counts("","IFIT1") # 404 not found
get_DE_counts("Bt.24795.1.A1_at","IFIT2")
plot_expression("Bt.24795.1.A1_at","IFIT2")
boxplot_expression("Bt.24795.1.A1_at","IFIT2")

get_DE_counts("Bt.4675.1.S1_a_at","MX1")
plot_expression("Bt.4675.1.S1_a_at","MX1")
boxplot_expression("Bt.4675.1.S1_a_at","MX1")

get_DE_counts("Bt.8143.1.S1_at","MX2")
plot_expression("Bt.8143.1.S1_at","MX2")
boxplot_expression("Bt.8143.1.S1_at","MX2")
##get_DE_counts("","IL27") # 404 not found

# NF-kB-dependent genes
plot_expression("Bt.12756.1.S1_at","TNF")
get_DE_counts("Bt.12756.1.S1_at","TNF")
boxplot_expression("Bt.12756.1.S1_at","TNF")

plot_expression("Bt.3686.1.S1_at","IL6")
get_DE_counts("Bt.3686.1.S1_at","IL6")
boxplot_expression("Bt.3686.1.S1_at","IL6")

plot_expression("Bt.9309.1.A1_at","NFKB1")
get_DE_counts("Bt.9309.1.A1_at","NFKB1") 
boxplot_expression("Bt.9309.1.A1_at","NFKB1")

plot_expression("Bt.20288.1.S1_at","NFKB2")
get_DE_counts("Bt.20288.1.S1_at","NFKB2")
boxplot_expression("Bt.20288.1.S1_at","NFKB2")


plot_expression("Bt.4388.1.S1_at","IL27RA")
plot_expression("Bt.26956.1.S1_at","NLRP1")
plot_expression("Bt.11586.1.S1_at","PYCARD")

###
# Table probe <-> gene.symbol #
###
# I want to set up a simple website with a form
# where one gives the gene.symbol
# and obtains the list of probes mapped to it
# and then can show the expression plots associated with each of those probes
setwd("C:\\Users\\krue\\Documents\\Kevin-Logs\\20130320_MDM\\Database")
tmp_probes = data.frame(ID=rownames(assayData(farms)$exprs))
tmp_probes$ID = rownames(assayData(farms)$exprs)
tmp_probes_geneSymbols = AnnotateTopTable(tmp_probes)[,c("ID","gene.symbols")]
write.table(x=tmp_probes_geneSymbols, file="genesymbols_probes.tab", quote=F, sep="\t", eol="\n", dec=".", row.names=F, col.names=F)
rm(tmp_probes, tmp_probes_geneSymbols)
# then imported in PhpMyAdmin

######################################
# Generate all boxplots of 5 animals #
######################################
setwd("C:\\Users\\krue\\Documents\\Kevin-Logs\\20130320_MDM\\Database\\boxplots")
# ALL probe sets !!
# Not just informative ones
# I plan to add the "informative" label in the database, to show alongside the boxplot of non-informative probe sets
tmp = data.frame(t(assayData(farms)$exprs))

tmp$group = NA #
tmp$group[grep("CN_0HR",rownames(tmp))] = "CN_0HR" #
tmp$group[grep("CN_2HR",rownames(tmp))] = "CN_2HR" #
tmp$group[grep("CN_6HR",rownames(tmp))] = "CN_6HR" #
tmp$group[grep("CN_24HR",rownames(tmp))] = "CN_24HR" #
tmp$group[grep("CN_25HR",rownames(tmp))] = "CN_24HR" #

tmp$group[grep("BCG_2HR",rownames(tmp))] = "BCG_2HR" #
tmp$group[grep("BCG_6HR",rownames(tmp))] = "BCG_6HR" #
tmp$group[grep("BCG_24HR",rownames(tmp))] = "BCG_24HR" #
tmp$group[grep("BCG_25HR",rownames(tmp))] = "BCG_24HR" #

tmp$group[grep("MAP_2HR",rownames(tmp))] = "MAP_2HR" #
tmp$group[grep("MAP_6HR",rownames(tmp))] = "MAP_6HR" #
tmp$group[grep("MAP_24HR",rownames(tmp))] = "MAP_24HR" #
tmp$group[grep("MAP_25HR",rownames(tmp))] = "MAP_24HR" #

tmp$group[grep("BOVIS_2HR",rownames(tmp))] = "BOVIS_2HR" #
tmp$group[grep("BOVIS_6HR",rownames(tmp))] = "BOVIS_6HR" #
tmp$group[grep("BOVIS_24HR",rownames(tmp))] = "BOVIS_24HR" #
tmp$group[grep("BOVIS_25HR",rownames(tmp))] = "BOVIS_24HR" #

tmp$group = as.factor(tmp$group) #
tmp$group = factor(tmp$group, groups) #

setwd("C:\\Users\\krue\\Documents\\Kevin-Logs\\20130320_MDM\\Database\\boxplots2")
for(tmpProbe in rownames(assayData(farms)$exprs))
{
  tmpGS = getSYMBOL(tmpProbe, data="bovine.db")
  png(filename=paste(tmpProbe,"-",tmpGS, ".png", sep=""), width=700, height=700)
  tryCatch({
    boxplot(tmp[,tmpProbe]~group, data=tmp, cex.axis=0.78, las=3, main=paste(tmpGS, "=", tmpProbe),
            col=c(rep("green",4), rep("blue",3), rep("brown",3), rep("red",3)))
  }, warning = function(w) {
    print(w)
  }, error = function(e) {
    boxplot(tmp[,gsub("-", "\\.", tmpProbe)]~group, data=tmp, cex.axis=0.78, las=3, main=paste(tmpGS, "=", tmpProbe),
            col=c(rep("green",4), rep("blue",3), rep("brown",3), rep("red",3)))
  }, finally = {
    dev.off()}
  )
}

###################
# Common DE genes #
###################
setwd("C:/Users/krue/Documents/Kevin-Logs/20130320_MDM/IPA-Overlap/")

# There are 415 and 541 DE gens common to all three infections at 2 hours and 6 hours, respectively
# The question is: how much do these overlap?
## Venn diagram of "common at 2hrs" and "common at 6 hours"


# first define if each informative probe is "common at 2hrs" 
common.df = data.frame(rownames=DE.count$ID)
common.df$ID = DE.count$ID
common.df$common2 = apply(DE.count[,c("MAP_2HR-CN_2HR", "BCG_2HR-CN_2HR", "BOVIS_2HR-CN_2HR")], MARGIN=1, FUN=all)
common.df$common6 = apply(DE.count[,c("MAP_6HR-CN_6HR", "BCG_6HR-CN_6HR", "BOVIS_6HR-CN_6HR")], MARGIN=1, FUN=all)
common.df$common24 = apply(DE.count[,c("MAP_24HR-CN_24HR", "BCG_24HR-CN_24HR", "BOVIS_24HR-CN_24HR")], MARGIN=1, FUN=all)

sum(common.df$common2) # 415
sum(common.df$common6) # 541
sum(common.df$common24) # 4

vennDiagram(common.df[,c("common2", "common6")])

# will export the genes common at 2 hpi for DAVID TFBS analysis ######################
setwd("C:/Users/krue/Documents/Kevin-Logs/20130320_MDM/IPA-Overlap/")
tmp.common = AnnotateTopTable(common.df)
tmp.common[tmp.common$common2,]
tmp.common[tmp.common$common2 & !isNA(tmp.common$gene.symbols),]
tmp.common[tmp.common$common2 & !isNA(tmp.common$gene.symbols),]$gene.symbols # 313 (redundant)
unique(tmp.common[tmp.common$common2 & !isNA(tmp.common$gene.symbols),]$gene.symbols) # 242 (unique)
write.table(unique(tmp.common[tmp.common$common2 & !isNA(tmp.common$gene.symbols),]$gene.symbols), file="genesymbols_common_2hpi.txt", quote=F, sep="\t", eol="\n", dec=".", row.names=F, col.names=F)
#common.df[common.df$common2,]$ID # 415
#write.table(common.df[common.df$common2,]$ID, file="probeIDs_common_2hpi.txt", quote=F, sep="\t", eol="\n", dec=".", row.names=F, col.names=F)
rm(tmp.common)

# NFKB and SP1 have enriched binding sites in the common_2hpi dataset
# It seems that SP1 is a subset of NFKB
# I want to check that
tmp.sP1 = read.csv(file="DAVID_SP1_TFBS_common_2hpi.txt", header=T, sep="\t", quote="")
tmp.NFKB = read.csv(file="DAVID_NFKB_TFBS_common_2hpi.txt", header=T, sep="\t", quote="")
tmp.union = data.frame(ID=union(tmp.sP1$ID, tmp.NFKB$ID))
tmp.union$SP1 = tmp.union$ID %in% tmp.sP1$ID
tmp.union$NFKB = tmp.union$ID %in% tmp.NFKB$ID
vennDiagram(object=tmp.union[,c("NFKB","SP1")])
rm(tmp.sP1, tmp.NFKB, tmp.union)

# same analysis for genes "common at 6hrs" ##############################
setwd("C:/Users/krue/Documents/Kevin-Logs/20130320_MDM/IPA-Overlap/")
unique(tmp.common[tmp.common$common6 & !isNA(tmp.common$gene.symbols),]$gene.symbols) # 328 unique gene symbols
write.table(unique(tmp.common[tmp.common$common6 & !isNA(tmp.common$gene.symbols),]$gene.symbols), file="genesymbols_common_6hpi.txt", quote=F, sep="\t", eol="\n", dec=".", row.names=F, col.names=F)

# In this case, there are just too many transcription factors which have enriched motifs in the uploaded dataset
# NFkB is still one of them.


############################################################
# Correlation of Diff Expr between mycobacterial treatments # 
############################################################
setwd("C:/Users/krue/Documents/Kevin-Logs/20130320_MDM/mycobacterial-correlation/")

# I want to plot the log2FC of mycobacteria-vs-control against one another
# To see how well the fold changes caused by one mycobacteria correlates with another
# Most of the boxplots show a similar direction of fold change for all mycobacteria
# although higher changes for BOVIS than others. The correlation should reflect this.
cn_cntrs = c("BOVIS_2HR-CN_2HR","BOVIS_6HR-CN_6HR","BOVIS_24HR-CN_24HR",
              "MAP_2HR-CN_2HR","MAP_6HR-CN_6HR","MAP_24HR-CN_24HR",
              "BCG_2HR-CN_2HR","BCG_6HR-CN_6HR","BCG_24HR-CN_24HR")
cn.FCs_df = data.frame(row.names=1:iProbes)
cn.FCs_df$ID=topTable(eb_filtered, coef= 1, number=iProbes, adjust.method="BH", sort.by="P")$ID
for(cntr in cn_cntrs){
  topT = topTable(eb_filtered, coef= cntr, number=iProbes, adjust.method="BH", sort.by="P")
  names(topT)[names(topT)=="logFC"] <- cntr
  cn.FCs_df = merge(x=cn.FCs_df, y=topT[c("ID",cntr)], by="ID", sort=T)
}
rm(topT)
write.table(cn.FCs_df, file = "logFC_vs-CN.txt", quote = FALSE, sep = "\t", eol = "\n", na = "NA", row.names = FALSE, col.names = TRUE)

# Correlation values ###
tmp.cor = cor(x=cn.FCs_df[,-1], method="pearson")
write.table(x=tmp.cor, file="pearson_correlation_matrix.txt", append=F, quote=F, sep="\t",
            eol="\n", dec=".", row.names=T, col.names=T)

# just out of curiosity the spearman one (non parametric)
tmp.cor = cor(x=cn.FCs_df[,-1], method="spearman")
write.table(x=tmp.cor, file="spearman_correlation_matrix.txt", append=F, quote=F, sep="\t",
            eol="\n", dec=".", row.names=T, col.names=T)
rm(tmp.cor)

# the pearson correlation is much better, which is a good sign. It means that there is an actual linear relationship between the 
# log2FC in the various conditions.
# On the other hand, the spearman returns lower correlation values, meaning that the order of probes based on their log2FC is not 
# as well conserved across treatments. This makes sense given the noisy pattern of logFC around the zero value.
# We are more interested in knowing whether a logFC value in one contrast correlates with the vale in other contrasts. There are so many
# probes with similar low logFC values that their order changes a lot between contrasts.

# Plotting ###
par(mfrow=c(1,3))
min(cn.FCs_df[,-1]) # -2.943241 is the min fold change
max(cn.FCs_df[,-1]) # 5.879359 is the max fold change (useful for xlim and ylim)
######################
# 2 hours ##
#####################
# BOV/MAP #
p <- qplot(cn.FCs_df$"BOVIS_2HR-CN_2HR", cn.FCs_df$"MAP_2HR-CN_2HR",
           xlab = "BOVIS_2HR-CN_2HR", ylab="MAP_2HR-CN_2HR",
           xlim=c(-3, 6), ylim=c(-3, 6))
# Shows the graph
p
# Plots the best fit line
p + stat_smooth(method="lm", se=TRUE)
coef(lm(`MAP_2HR-CN_2HR`~`BOVIS_2HR-CN_2HR`, data = cn.FCs_df))
#  (Intercept) `BOVIS_2HR-CN_2HR` 
#  0.003202317        0.526880374 
######################
# BOV/BCG #
p <- qplot(cn.FCs_df$"BOVIS_2HR-CN_2HR", cn.FCs_df$"BCG_2HR-CN_2HR",
           xlab = "BOVIS_2HR-CN_2HR", ylab="BCG_2HR-CN_2HR",
           xlim=c(-3, 6), ylim=c(-3, 6))
# Shows the graph
p
# Plots the best fit line
p + stat_smooth(method="lm", se=TRUE)
coef(lm(`BCG_2HR-CN_2HR`~`BOVIS_2HR-CN_2HR`, data = cn.FCs_df))
# (Intercept) `BOVIS_2HR-CN_2HR` 
# 0.0004648545       0.5421746591
########################
# MAP/BCG # 
p <- qplot(cn.FCs_df$"MAP_2HR-CN_2HR", cn.FCs_df$"BCG_2HR-CN_2HR",
           xlab = "MAP_2HR-CN_2HR", ylab="BCG_2HR-CN_2HR",
           xlim=c(-3, 6), ylim=c(-3, 6))

# Shows the graph
p
# Plots the best fit line
p + stat_smooth(method="lm", se=TRUE)
coef(lm(`BCG_2HR-CN_2HR`~`MAP_2HR-CN_2HR`, data = cn.FCs_df))
# (Intercept) `MAP_2HR-CN_2HR` 
# 0.003664996      1.009780540

########################
# 6 hours
########################
# BOV/MAP #
p <- qplot(cn.FCs_df$"BOVIS_6HR-CN_6HR", cn.FCs_df$"MAP_6HR-CN_6HR",
           xlab = "BOVIS_6HR-CN_6HR", ylab="MAP_6HR-CN_6HR",
           xlim=c(-3, 6), ylim=c(-3, 6))
# Shows the graph
p
# Plots the best fit line
p + stat_smooth(method="lm", se=TRUE)
coef(lm(`MAP_6HR-CN_6HR`~`BOVIS_6HR-CN_6HR`, data = cn.FCs_df))
#  (Intercept) `BOVIS_2HR-CN_2HR` 
#  0.003202317        0.526880374 
######################
# BOV/BCG #
p <- qplot(cn.FCs_df$"BOVIS_6HR-CN_6HR", cn.FCs_df$"BCG_6HR-CN_6HR",
           xlab = "BOVIS_6HR-CN_6HR", ylab="BCG_6HR-CN_6HR",
           xlim=c(-3, 6), ylim=c(-3, 6))
# Shows the graph
p
# Plots the best fit line
p + stat_smooth(method="lm", se=TRUE)
coef(lm(`BCG_6HR-CN_6HR`~`BOVIS_6HR-CN_6HR`, data = cn.FCs_df))
# (Intercept) `BOVIS_6HR-CN_6HR` 
# 0.0004648545       0.5421746591
########################
# MAP/BCG # 
p <- qplot(cn.FCs_df$"MAP_6HR-CN_6HR", cn.FCs_df$"BCG_6HR-CN_6HR",
           xlab = "MAP_6HR-CN_6HR", ylab="BCG_6HR-CN_6HR",
           xlim=c(-3, 6), ylim=c(-3, 6))

# Shows the graph
p
# Plots the best fit line
p + stat_smooth(method="lm", se=TRUE)
coef(lm(`BCG_6HR-CN_6HR`~`MAP_6HR-CN_6HR`, data = cn.FCs_df))
# (Intercept) `MAP_6HR-CN_6HR` 
# 0.003664996      1.009780540

########################
# 24 hours
########################
# BOV/MAP #
p <- qplot(cn.FCs_df$"BOVIS_24HR-CN_24HR", cn.FCs_df$"MAP_24HR-CN_24HR",
           xlab = "BOVIS_24HR-CN_24HR", ylab="MAP_24HR-CN_24HR",
           xlim=c(-3, 6), ylim=c(-3, 6))
# Shows the graph
p
# Plots the best fit line
p + stat_smooth(method="lm", se=TRUE)
coef(lm(`MAP_24HR-CN_24HR`~`BOVIS_24HR-CN_24HR`, data = cn.FCs_df))
#  (Intercept) `BOVIS_2HR-CN_2HR` 
#  0.003202317        0.526880374 
######################
# BOV/BCG #
p <- qplot(cn.FCs_df$"BOVIS_24HR-CN_24HR", cn.FCs_df$"BCG_24HR-CN_24HR",
           xlab = "BOVIS_24HR-CN_24HR", ylab="BCG_24HR-CN_24HR",
           xlim=c(-3, 6), ylim=c(-3, 6))
# Shows the graph
p
# Plots the best fit line
p + stat_smooth(method="lm", se=TRUE)
coef(lm(`BCG_24HR-CN_24HR`~`BOVIS_24HR-CN_24HR`, data = cn.FCs_df))
# (Intercept) `BOVIS_6HR-CN_6HR` 
# 0.0004648545       0.5421746591
########################
# MAP/BCG # 
p <- qplot(cn.FCs_df$"MAP_24HR-CN_24HR", cn.FCs_df$"BCG_24HR-CN_24HR",
           xlab = "MAP_24HR-CN_24HR", ylab="BCG_24HR-CN_24HR",
           xlim=c(-3, 6), ylim=c(-3, 6))

# Shows the graph
p
# Plots the best fit line
p + stat_smooth(method="lm", se=TRUE)
coef(lm(`BCG_24HR-CN_24HR`~`MAP_24HR-CN_24HR`, data = cn.FCs_df))
# (Intercept) `MAP_6HR-CN_6HR` 
# 0.003664996      1.009780540

rm(p)
##########################
# Filter on fold changes #
##########################
cn.FCs_df[cn.FCs_df$"MAP_24HR-CN_24HR" > 2,]
cn.FCs_df[cn.FCs_df$"BCG_24HR-CN_24HR" > 2,]


#######################
# Test of correlation #
#######################
# How come all the correlation are stated significant, even with R as low as 0.45?
test = cor.test(cn.FCs_df$`MAP_24HR-CN_24HR`,cn.FCs_df$`BCG_2HR-CN_2HR`,)
test
rm(test)
# The significance call tells that the correlation coefficient is signif different from zero
# (even weak, there is a relationship between the two values)
#
# Note that stronger correlation are found between fold changes at the same time point
# obviously, those samples "correlate"


#############
# Gene shared down-regulated (2hrs)
#############

setwd("C:/Users/krue/Documents/Kevin-Logs/20130320_MDM/based-on-Liz-results")

##### Extract from code for macrosys5 database 
# initialises the data frame with the appropriate number of rows
tmp_probes = data.frame(row.names=1:iProbes) 

# gets the IDs of the informative probe sets
tmp_probes$ID = topTable(eb_filtered, coef=1, number=iProbes, adjust.method="BH")$ID

# Progressively merges the ID table with logFC and adp.p.value information from each contrast
for(cntr in cntrts){
  topT = topTable(eb_filtered, coef= cntr, number=iProbes, adjust.method="BH")[,c("ID","logFC", "adj.P.Val")]
  colnames(topT)[colnames(topT)=="logFC"] <- paste("logFC.", cntr, sep="")
  colnames(topT)[colnames(topT)=="adj.P.Val"] <- paste("adj.P.Val.", cntr, sep="")
  tmp_probes = merge(x=tmp_probes, y=topT, by="ID", sort=TRUE)
  rm(topT)
}

##### Need annotations
tmp.annotFC = AnnotateTopTable(df=tmp_probes)

####### Filter for probe shared down-regulated
tmp = tmp.annotFC[tmp.annotFC["logFC.BOVIS_2HR-CN_2HR"] < 0 &
              tmp.annotFC["logFC.MAP_2HR-CN_2HR"] < 0 &
              tmp.annotFC["logFC.BCG_2HR-CN_2HR"] < 0 &
              tmp.annotFC["adj.P.Val.BOVIS_2HR-CN_2HR"] < 0.05 &
              tmp.annotFC["adj.P.Val.BCG_2HR-CN_2HR"] < 0.05 &
              tmp.annotFC["adj.P.Val.MAP_2HR-CN_2HR"] < 0.05 , ]$ID
tmp.filtered = data.frame(ID=tmp)
tmp.filtered$ID = tmp
tmp.filtered = AnnotateTopTable(tmp.filtered)
write.table(x=tmp.filtered, file="shared_2hpi-downregulated.txt", append=F, quote=F, sep="\t",
            eol="\n", dec=".", row.names=F, col.names=T)
tmp.downs = tmp.filtered

####### Filter for probe shared up-regulated
tmp = tmp.annotFC[tmp.annotFC["logFC.BOVIS_2HR-CN_2HR"] > 0 &
                    tmp.annotFC["logFC.MAP_2HR-CN_2HR"] > 0 &
                    tmp.annotFC["logFC.BCG_2HR-CN_2HR"] > 0 &
                    tmp.annotFC["adj.P.Val.BOVIS_2HR-CN_2HR"] < 0.05 &
                    tmp.annotFC["adj.P.Val.BCG_2HR-CN_2HR"] < 0.05 &
                    tmp.annotFC["adj.P.Val.MAP_2HR-CN_2HR"] < 0.05 , ]$ID
tmp.filtered = data.frame(ID=tmp)
tmp.filtered$ID = tmp
tmp.filtered = AnnotateTopTable(tmp.filtered)
write.table(x=tmp.filtered, file="shared_2hpi-upregulated.txt", append=F, quote=F, sep="\t",
            eol="\n", dec=".", row.names=F, col.names=T)

tmp.ups = tmp.filtered
####### Merge probes shared up- and down- regulated
tmp.downs$UpDown = -1
tmp.ups$UpDown = 1
tmp.merged = rbind(tmp.downs, tmp.ups)
write.table(x=tmp.merged, file="shared_2hpi-DE.txt", append=F, quote=F, sep="\t",
            eol="\n", dec=".", row.names=F, col.names=T)





#############
# Gene shared down-regulated (6hrs)
#############

setwd("C:/Users/krue/Documents/Kevin-Logs/20130320_MDM/based-on-Liz-results")

# use tmp.annotFC containing the logFC and p-values

####### Filter for probe shared down-regulated
tmp = tmp.annotFC[tmp.annotFC["logFC.BOVIS_6HR-CN_6HR"] < 0 &
                    tmp.annotFC["logFC.MAP_6HR-CN_6HR"] < 0 &
                    tmp.annotFC["logFC.BCG_6HR-CN_6HR"] < 0 &
                    tmp.annotFC["adj.P.Val.BOVIS_6HR-CN_6HR"] < 0.05 &
                    tmp.annotFC["adj.P.Val.BCG_6HR-CN_6HR"] < 0.05 &
                    tmp.annotFC["adj.P.Val.MAP_6HR-CN_6HR"] < 0.05 , ]$ID
tmp.filtered = data.frame(ID=tmp)
tmp.filtered$ID = tmp # 123 down
tmp.filtered = AnnotateTopTable(tmp.filtered)
write.table(x=tmp.filtered, file="shared_6hpi-downregulated.txt", append=F, quote=F, sep="\t",
            eol="\n", dec=".", row.names=F, col.names=T)
tmp.downs = tmp.filtered

####### Filter for probe shared up-regulated
tmp = tmp.annotFC[tmp.annotFC["logFC.BOVIS_6HR-CN_6HR"] > 0 &
                    tmp.annotFC["logFC.MAP_6HR-CN_6HR"] > 0 &
                    tmp.annotFC["logFC.BCG_6HR-CN_6HR"] > 0 &
                    tmp.annotFC["adj.P.Val.BOVIS_6HR-CN_6HR"] < 0.05 &
                    tmp.annotFC["adj.P.Val.BCG_6HR-CN_6HR"] < 0.05 &
                    tmp.annotFC["adj.P.Val.MAP_6HR-CN_6HR"] < 0.05 , ]$ID
tmp.filtered = data.frame(ID=tmp)
tmp.filtered$ID = tmp # 421 up
tmp.filtered = AnnotateTopTable(tmp.filtered)
write.table(x=tmp.filtered, file="shared_6hpi-upregulated.txt", append=F, quote=F, sep="\t",
            eol="\n", dec=".", row.names=F, col.names=T)

tmp.ups = tmp.filtered
####### Merge probes shared up- and down- regulated
tmp.downs$UpDown = -1
tmp.ups$UpDown = 1
tmp.merged = rbind(tmp.downs, tmp.ups)
write.table(x=tmp.merged, file="shared_6hpi-DE.txt", append=F, quote=F, sep="\t",
            eol="\n", dec=".", row.names=F, col.names=T)

#############
# Heatmap + dendrograms
#############
heatplot(dataset=farms_informative,
         dualScale=FALSE, scale="none",
         cexCol=0.5, cexRow=0.5)

heatplot(dataset=farms_informative[1:100,], dualScale=FALSE, # faster for a subset of genes.
         scale="none", cexCol=0.5)
# note that both "Affx" and "Bt" probes are present among the informative ones
# We would probably need to filter for:
# probes with large variance (see Teles 2013)
# probes filtered out because ANOVA only returns animal as a significant variance factor.

#############
# ANOVA(s) (The code below does not do an ANOVA. I tried to get there but couln't)
#############
#tmp.aov = aov(formula=)
#head(assayData(farms_informative)$exprs)
tmp.design = createDesignMatrix(farms_informative)
tmp.contM = createContrastMatrix(farms_informative)

# plotMDS.default(x=farms)
plotMDS.default(x=farms_informative)  # exact same plot as previous line

# Lines below copied from Paul's help
# for DE genes and interaction effects

# prepare the data for Puma
farms_informative2 = farms_informative  # duplicate the data because the puma functions are more restrictive about the data layout
farms_informative2 = farms_informative2[,which(farms_informative2$TimePoint != "0HR")] # the zero time point prevent comparison across treatment and time points
pData(farms_informative2) = pData(farms_informative2)[pData(farms_informative2)$TimePoint != "0HR",2:4] # the extra columns are confusing for eBayes which does not need them
# the actual analysis
design2 = createDesignMatrix(farms_informative2)
fit2 = lmFit(farms_informative2, design2)
contM2 = createContrastMatrix(farms_informative2)
fitC2 = contrasts.fit(fit2, contM2)
#fitC2 = eBayes(fitC2)


###########
# Heatplot like Teles 2013 (still cluster sby animal and time ! Instead of treatment and time)
############
setwd("C:/Users/krue/Documents/Kevin-Logs/20130320_MDM/based-on-teles 2013")

tmp = as.data.frame(assayData(farms_informative)$exprs) # expr values matrix (row: probe, col: sample)
tmp$var = apply(tmp, MARGIN=1, FUN=var) # calculates the variance for each probe between all samples

sum(tmp$var >= 1) # how many probes have variance >= 1? 261
which(tmp$var >= 1) # which rows indices correspond to the above probes?
tmp.var1plus = rownames(tmp[which(tmp$var >= 1),]) # which probes IDs are those?
tmp.farmsInfo_var1plus = farms_informative[tmp.var1plus,] # expression data for those probes
tmp.farmsInfo_var1plus # how many of them are informative? 261

heatplot(dataset=tmp.farmsInfo_var1plus, # Heatmap and dendrogram based on those genes
         dualScale=FALSE, scale="none",
         cexCol=0.5, cexRow=0.5)

tmp.df = data.frame(ID=tmp.var1plus) # annotate those probes
tmp.df$ID = tmp.var1plus
tmp.df = AnnotateTopTable(tmp.df)

write.table(x=tmp.df, file="261genes_informative_var1plus.txt", append=F, quote=F, sep="\t",
            eol="\n", dec=".", row.names=F, col.names=T)



################
# PUMA package #
################
## library(PUMA)
#
# With Paul's Advice
# Given that the samples keep clustering by animal rather than treatment
# He suggested to call DE genes blocking on a different factor than animal to see th effect of animal on DE call
# There may be genes called DE between animals, forcing them to cluster apart in the dendrogram
#
# ##### Animal + Timepoint #########
farms2 = farms  # duplicate the data because the puma functions are more restrictive about the data layout
farms2 = farms2[,which(farms2$TimePoint != "0HR")] # the zero time point prevent comparison across treatment and time points
pData(farms2) = pData(farms)[pData(farms)$TimePoint != "0HR",c(2,4)] # keeps animal and timepoint as explanatory factors
# the actual analysis
design = createDesignMatrix(farms2)
fit = lmFit(farms2, design)
contM = createContrastMatrix(farms2)
fitC = contrasts.fit(fit,contM)
fitC = eBayes(fitC)
# Get the number DE calls at q-value 0.05
esClas = decideTests(fitC,p.value=0.05,adjust.method="BH")
# Get the total count of DE genes for each contrast
colSums(abs(esClas))
#write.table(x=colSums(abs(esClas)), file="PUMA-Animal-Timepoint-DEcalls.txt", append=FALSE, quote=FALSE, sep="\t", eol="\n", dec=".", col.names="DE.counts")


# ##### Animal + Treatment #########
farms2 = farms  # duplicate the data because the puma functions are more restrictive about the data layout
farms2 = farms2[,which(farms2$TimePoint != "0HR")] # the zero time point prevent comparison across treatment and time points
pData(farms2) = pData(farms)[pData(farms)$TimePoint != "0HR", 2:3] # keeps animal and treatment as explanatory factors
# the actual analysis
design = createDesignMatrix(farms2)
fit = lmFit(farms2, design)
contM = createContrastMatrix(farms2)
fitC = contrasts.fit(fit,contM)
fitC = eBayes(fitC)
# Get the number DE calls at q-value 0.05
esClas = decideTests(fitC,p.value=0.05,adjust.method="BH")
# Get the total count of DE genes for each contrast
colSums(abs(esClas))
#write.table(x=colSums(abs(esClas)), file="PUMA-Animal-Timepoint-DEcalls.txt", append=FALSE, quote=FALSE, sep="\t", eol="\n", dec=".", col.names="DE.counts")

# Get genes differentially expressed between two treatments within one animal (merging all time points)
tmp.names = names(esClas[esClas[,"721.MAP_vs_721.BOVIS"] != 0,][,"721.MAP_vs_721.BOVIS"]) # get the probeIDs of 94 DE genes of animal 721 MAP-BOVIS
boxplot(exprs(farms_informative2)[tmp.names[1],]~pData(farms_informative2)$Animal+pData(farms_informative2)$Treatment) # boxplot of expression values for first of those probes
tmp.df = data.frame(ID=tmp.names)# annotate the above probeIDs
tmp.df$ID = tmp.names
tmp.df = AnnotateTopTable(tmp.df)

write.table(x=tmp.df, file="94genes_informative_DE_A-721_BOVIS-MAP.txt", append=F, quote=F, sep="\t",
            eol="\n", dec=".", row.names=F, col.names=T)

# In the first 5 annotated genes, RAB7B can be found
# very interestingly: http://www.genecards.org/cgi-bin/carddisp.pl?gene=RAB7B
boxplot(exprs(farms_informative2)["Bt.11413.1.S1_at",]~pData(farms_informative2)$Animal+pData(farms_informative2)$Treatment,
        las= 2, cex.axis=0.7, main="RAB7B = Bt.11413.1.S1_at") # boxplot of expression values for first of those probes


###### Strong animal-influenced genes ###
# Paul: Pick animals 706 and 721 (more or less randomyl picked) 
# Filter for genes DE in "Animal_706_vs_721" column (+/- 1) in Animal+Timepoint and Animal+Treatment analysis
# Boxplot some of the highest p-values genes (esClas of either analysis)
#

# ##### Animal + Timepoint #########
farms2 = farms  # duplicate the data because the puma functions are more restrictive about the data layout
farms2 = farms2[,which(farms2$TimePoint != "0HR")] # the zero time point prevent comparison across treatment and time points
pData(farms2) = pData(farms)[pData(farms)$TimePoint != "0HR",c(2,4)] # keeps animal and timepoint as explanatory factors
# the actual analysis
design = createDesignMatrix(farms2)
fit = lmFit(farms2, design)
contM = createContrastMatrix(farms2)
fitC1 = contrasts.fit(fit,contM)
fitC1 = eBayes(fitC1)
# Get the number DE calls at q-value 0.05
esClas1 = decideTests(fitC1,p.value=0.05,adjust.method="BH")

tmp.names1 = names(esClas1[esClas1[,"Animal_706_vs_721"] != 0,][,"Animal_706_vs_721"]) ### Filter for genes DE in "Animal_706_vs_721" column (+/- 1) in Animal+Timepoint

# ##### Animal + Treatment #########
farms2 = farms  # duplicate the data because the puma functions are more restrictive about the data layout
farms2 = farms2[,which(farms2$TimePoint != "0HR")] # the zero time point prevent comparison across treatment and time points
pData(farms2) = pData(farms)[pData(farms)$TimePoint != "0HR", 2:3] # keeps animal and treatment as explanatory factors
# the actual analysis
design = createDesignMatrix(farms2)
fit = lmFit(farms2, design)
contM = createContrastMatrix(farms2)
fitC2 = contrasts.fit(fit,contM)
fitC2 = eBayes(fitC2)
# Get the number DE calls at q-value 0.05
esClas2 = decideTests(fitC2,p.value=0.05,adjust.method="BH")

tmp.names2 = names(esClas2[esClas2[,"Animal_706_vs_721"] != 0,][,"Animal_706_vs_721"]) ### Filter for genes DE in "Animal_706_vs_721" column (+/- 1) in Animal+Timepoint

tmp.names.inter = intersect(tmp.names1, tmp.names2) # Probes DE between Animal 706 and 721, both obtained from Time and Treatment blockings
tmp.df = data.frame(ID=tmp.names.inter) # turn it into a dataset
tmp.df$ID = tmp.names.inter # make it easy to annotate

# obtain the  p-value for each probe for each comparison
fitC1.df = data.frame(fitC1)  # prepare Animal+Time data as frame
fitC1.df$ID = rownames(fitC1.df) # make it easy to merge
tmp.df = merge(x=tmp.df, y=fitC1.df[,c("p.value.Animal_706_vs_721","ID")], by.x="ID", by.y="ID") # add the column to the dateset
names(tmp.df)[names(tmp.df) == "p.value.Animal_706_vs_721"] = "Animal.Time.p.value.Animal_706_vs_721" # rename the column because the next will have the same name in Anmal+Treatment

fitC2.df = data.frame(fitC2)  # prepare Animal+Treatment data as frame
fitC2.df$ID = rownames(fitC2.df) # make it easy to merge
tmp.df = merge(x=tmp.df, y=fitC2.df[,c("p.value.Animal_706_vs_721","ID")], by.x="ID", by.y="ID") # add the column to the dateset
names(tmp.df)[names(tmp.df) == "p.value.Animal_706_vs_721"] = "Animal.Treatment.p.value.Animal_706_vs_721" # rename the column because the next will have the same name in Anmal+Treatment

# average log10 of both p-values for ranking purpose
tmp.df$meanlogPvalue = (log10(tmp.df$"Animal.Time.p.value.Animal_706_vs_721")+log10(tmp.df$Animal.Treatment.p.value.Animal_706_vs_721))/2
tmp.df = tmp.df[order(tmp.df$meanlogPvalue),] # ranking

tmp.df = AnnotateTopTable(tmp.df)

#### At last: the boxplots !!!
# The first annotated DE gene is a BOLA (MHC-II)
# ... just like the next 10
# They all have an average log p value between -30 and -60
boxplot(exprs(farms_informative2)["Bt.4751.1.S1_a_at",]~pData(farms_informative2)$Treatment+pData(farms_informative2)$Animal,
        las= 2, cex.axis=0.7, main="BOLA-DQA2 = Bt.4751.1.S1_a_at") # boxplot of expression values for first of those probes
boxplot(exprs(farms_informative2)["Bt.4751.1.S1_a_at",]~pData(farms_informative2)$TimePoint+pData(farms_informative2)$Animal, # Beware these lines may not work depending on the factors kept in the pData ~30-40 lines above
        las= 2, cex.axis=0.7, main="BOLA-DQA2 = Bt.4751.1.S1_a_at") # boxplot of expression values for first of those probes


write.table(x=tmp.df, file="2171probes_informative_DE_A706-721.txt", append=F, quote=F, sep="\t",
            eol="\n", dec=".", row.names=F, col.names=T)

# Just checking a non-BOLA gene
boxplot(exprs(farms_informative2)["Bt.23296.1.S1_at",]~pData(farms_informative2)$Treatment+pData(farms_informative2)$Animal,
        las= 2, cex.axis=0.7, main="ATP2B1 = Bt.23296.1.S1_at") # boxplot of expression values for first of those probes

#IFNG has animal effect too
boxplot(exprs(farms_informative2)["Bt.188.1.S1_at",]~pData(farms_informative2)$Treatment+pData(farms_informative2)$Animal,
        las= 2, cex.axis=0.7, main="IFNG = Bt.188.1.S1_at") # boxplot of expression values for first of those probes



#
#
#
#
which(esClas2[,"Int__Animal_724R.727R_vs_Treatment_BCG.CN"] != 0) ### Filter for genes DE in "Int__Animal_724R.727R_vs_Treatment_BCG.CN" column (+/- 1) in Animal+Timepoint

boxplot(exprs(farms_informative2)["Bt.12825.1.S1_at",]~pData(farms_informative2)$Treatment+pData(farms_informative2)$Animal,
        las= 2, cex.axis=0.7, main="ACTA2 = Bt.12825.1.S1_at") # boxplot of expression values for first of those probes


# Showing the data to David Magee (can be useful later for me)

#head(fitC1$coefficients)
#fitC2$coefficients
#head(fitC1$coefficients[,"Animal_706_vs_721"])
#head(fitC1$p.value[,"Animal_706_vs_721"])
#fitC1$coefficients["Bt.4751.1.S1_a_at","Animal_706_vs_721"]


#############
# Removing the animal effect 
#############

# I want to identify the most likely value representing the animal effect
# which needs to be removed from the normalised intensity to get rid of that effect
# NOTE: there is an animal effect measured when dropping the time factor
#       and another when dropping the treatment factor.
#       I will try one at a time, before trying to combine them both which is the ultimate goal
# For this, I will store the distance from the mean for each probe for each sample.
# Then I will correlated different likely relevant metrics from fitC with the above distance.
# Hopefully I will find the metric which represents teh animal factor best
# Then I might have to transform this value (sqrt, squared, ...) to obtain a meaningful value to substract from the normalised intensity. 

exprs2 = assayData(farms2)$exprs # extract the normalised intensity for each probe for each sample
means = matrix(data=apply(exprs2, MARGIN=1, mean), nrow=length(rownames(exprs2)), dimnames=list(rownames(exprs2), "mean"))  # mean expression for each gene 
distances2 = sweep(x=exprs2, MARGIN=1, STATS=means, FUN="-") # distance from the mean for each sample for each gene

# checking the distance for top BOLA (Bt.4751.1.S1_a_at)
distances2["Bt.4751.1.S1_a_at",]
abline(h=means["Bt.4751.1.S1_a_at",])

# plots of metrics vs. above distance from the mean
head(fitC1$coefficients) # let'#s try Animal+Time dataset first
colnames(fitC1$coefficients)
fitC1$coefficients["Bt.4751.1.S1_a_at",]
plot(fitC1$coefficients~distances2)

typeof(fitC1) # says "list", but is MArrayLM, see http://svitsrv25.epfl.ch/R-doc/library/limma/html/marraylm.html
head(fitC1$genes)

fitC1$proportion

head(fitC1$stdev.unscaled)
fitC1$stdev.unscaled["Bt.4751.1.S1_a_at",]

# In fact, given the colnames found for coefficient, and other data structures, I don't think the can find a value measuring the animal effect 
# per animal. Indeed, we need to measure a signle value per animal and per probe, representing the effect of that particular animal for 
# for that particular probe. And then we could subtract that value from the normalised intensity of that probe in that animal.
# To this end, I believe we have to run an ANOVA (Animal+Time+Treatment) for each probe, and there we should find  

#########
# Per probe ANOVA strategy
########
#  I thought that if I run an ANOVA for each gene, griouping values by giving animal, treatment and time information
# then the ANOVA would spit out the participation of each animal.
# I extracted the intensities for the top BOLA gene (archetype of animal effect)
# And gave the infromation about the different factors.
# It did confirm that thet animal effect outperforms all others, but absence of replicates prevented p-values to be calculated
# Moreover this only computes the variance between and within groups, but not the per-animal effect coefficient. 

# In short, the ANOVA approach confirms what is found by PUMA, without the estimate of factor and interactions effects
# See below

# Extracting the data for the top BOLA gene (Bt.4751.1.S1_a_at)

test.measures = assayData(farms_informative2)$exprs["Bt.4751.1.S1_a_at",]
test.measures

test.animal = as.factor(rep(c("706","716","721","724R","727R"),12))
test.animal
test.treatment = as.factor(c(rep("CN",15),rep("BCG",15),rep("MAP",15),rep("BOVIS",15)))
test.treatment
test.time = as.factor(rep(c(rep("2HR",5),rep("6HR",5),rep("24HR",5)),4))
test.time
cbind(test.measures, test.animal, test.time, test.treatment)


fit = lm(test.measures ~ test.animal*test.treatment*test.time)
aov = anova(fit)
aov

# What we are missing to estimate EACH animal effect is in fact another layer of replicates for fixed animal+time+treatment
# Paul suggested that I clone the normalised intensities, and inject it in PUMA as if we had those replicates.
# This is expected to absolutely screw the p-values estimated, but the coefficients of the different effects shouldn't change from the values
# observed when dropping a factor
#
farms3 = farms
pData(farms3)
str(farms3)
head(assayData(farms3)$exprs)


#Paul's pipeline
#Creating a chimeric dataset with fake replicates
eset<-exprs(farms3) # store original intensities in eset
eset2<-eset # clones a second expr set
colnames(eset2)<-paste(colnames(eset2),"1",sep=".") # renames the colnames as replicates
eset3<-eset # clones a third expr set
colnames(eset3)<-paste(colnames(eset3),"2",sep=".") # renames the colnames as replicates
esetAll<-cbind(eset,eset2,eset3) # mergees the three expr sets in a single one including the three cloned dataset
rm(eset,eset2,eset3)

#Creating a chimeric phenodata structure accompanying the above dataset
pdata<-pData(farms3) # store original phenodata in pdata
pdata2<-pdata # creates a clone pdata
rownames(pdata2)<-paste(rownames(pdata),"1",sep=".") # renames the colnames as replicates
pdata3<-pdata # creates a clone pdata
rownames(pdata3)<-paste(rownames(pdata),"2",sep="." ) # renames the colnames as replicates
pdataAll<-rbind(pdata,pdata2,pdata3) # mergees the threephenodata in a single one including the three cloned dataset
rm(pdata,pdata2,pdata3)

# Create a chimeric ExpressionSet data structure binding the above two structures together
chimeric.pData<-new("AnnotatedDataFrame", data=pdataAll)
chimeric.eSet<-new("ExpressionSet",exprs=esetAll, phenoData=chimeric.pData)

# Hang on, the time zero samples are present in those, we need to get rid of them
chimeric.eSet = chimeric.eSet[,which(chimeric.eSet$TimePoint != "0HR")] # the zero time point prevent comparison across treatment and time points
pData(chimeric.eSet) = pData(chimeric.eSet)[pData(chimeric.eSet)$TimePoint != "0HR", 2:4] # keeps animal and treatment as explanatory factors



# PUMA analysis
chimeric.design = createDesignMatrix(chimeric.eSet)
chimeric.fit = lmFit(chimeric.eSet, chimeric.design)
chimeric.contM = createContrastMatrix(chimeric.eSet)
chimeric.fitC = contrasts.fit(chimeric.fit,chimeric.contM)
chimeric.fitC = eBayes(chimeric.fitC)
# Get the number DE calls at q-value 0.05
chimeric.esClas = decideTests(chimeric.fitC,p.value=0.05,adjust.method="BH")

colnames(chimeric.esClas)
colSums(abs(esClas))

# Alright, it didn't work but I know what's wrong
# Again, we got an estimates from contrasts between the three factors
# but this approach will NEVER give us an estimate of EACH animal's impact on the measures
#
# The most intuititive approach which will UNDOUBTEDLY give us the expected result (and also the most reasonable approach)
# is to to compute for each animal the mean distance of each probe measurement from the mean of all samples for that probe
# This will give the deviation of each animal from the population mean, which we can then subtract/add to bring those animals back to a similar level
# without affecting the effect of time and treatments.

# For the publication interpretation, Note that those genes with an animal effect may influence the rest of the genes by their presence/absence in the sample 
# This list of animal-specific gene patterns may be just as valuable as the list of genes driven by treatment and time.
# I am subtracting the animal effect to identify more accurately the list of gene affected by time and treatment.
# But it will be just as useful to investigate those genes largely affected by the animal factor.

###########
# THE PROPER CORRECTION FOR ANIMAL INFLUENCE
###########

# Just before I start, Paul suggested I have a gene with low animal effect to compare before/after normalisation
# I picked the least significant (but still p<0.05) gene with an animal effect (use it after the proper normalisation has worked)
boxplot(exprs(farms_informative2)["Bt.25051.1.A1_at",]~pData(farms_informative2)$Treatment+pData(farms_informative2)$Animal,
        las= 2, cex.axis=0.7, main="TTC14 = Bt.25051.1.A1_at") # boxplot of expression values for first of those probes
abline(v=8.5)
abline(v=12.5)
abline(h=means["Bt.25051.1.A1_at",])


# FYI the dataset used below is farms2 generated as such
farms2 = farms
farms2 = farms2[,which(farms2$TimePoint != "0HR")] # the zero time point prevent comparison across treatment and time points
pData(farms2) = pData(farms)[pData(farms)$TimePoint != "0HR", 2:3] # keeps animal and treatment as explanatory factors
# Now
# 0.1 = calculate mean for each probe in the farms2 dataset
#means = matrix(data=apply(assayData(farms2)$exprs, MARGIN=1, mean), nrow=length(rownames(exprs2)), dimnames=list(rownames(exprs2), "mean"))  # mean expression for each gene 
means = rowMeans(assayData(farms2)$exprs)
# 0.2 = calculate mean for each animal for each probe 
animals = rep(c("706","716","721","724R","727R"),12)
#tapply(t(assayData(farms2)$exprs[1,]),INDEX=pData(farms2)$Animal, FUN="mean")
f1<-function(x){
  return(tapply(x,INDEX=pData(farms2)$Animal, FUN="mean"))
}
animals.means = t(apply(assayData(farms2)$exprs, 1, f1))
# 0.1+0.2 -> 1 = calculate distance for each animal for each probe from the overall mean value for that probe
f2 = function(x)
{
  return(x - means)
}
animals.distances = apply(X=animals.means, MARGIN=2, FUN="f2") # woohoooo, this is the animal influence !!!
animals.distances["Bt.25051.1.A1_at",] # TTC14
animals.distances["Bt.4751.1.S1_a_at",] #  top BOLA
head(animals.distances)

# 1+farms2 -> 2 = correct farms2 into farms2.animal.corrected by subtracting the distance for each animal for each probe
# I will build a matrix (24128*5) copying the distance of the animal from the overall mean for each gene depending on the sample column
# From left to right in the exprs matrix, if the sample is 706_BOVIS_2HR, then I cbind the 706 distance, etc
tmp.subtract = matrix(nrow=24128, dimnames=list(rownames(animals.distances)))
for(animal in pData(farms2)$Animal)
{
  tmp.subtract = cbind(tmp.subtract, animals.distances[,animal])
}
tmp.subtract = tmp.subtract[,-1]
colnames(tmp.subtract) = rownames(pData(farms2)) # tadaaaaaaaaa, here is the matrix to subtract from the normalised expression!

tmp.subtract["Bt.4751.1.S1_a_at",] #  top BOLA
#After that, I will simply subtract this matrix to the exprs one

# create a new expression set with normalised intensities CORRECTED for the animal effect
tmp.eset = (assayData(farms2)$exprs - tmp.subtract)

pdata<-new("AnnotatedDataFrame", data=pData(farms)[pData(farms)$TimePoint != "0HR",])
farms2.animal.corrected <- new("ExpressionSet",exprs=tmp.eset,phenoData=pdata, se.exprs=assayData(farms2)$se.exprs)

# Comparison of BOLA before/after
animals.abline = function()
{
  for(i in 0:3){
    abline(v=4*i+4.5)
  }
}
par(mfrow=c(2,2))
boxplot(exprs(farms_informative2)["Bt.4751.1.S1_a_at",]~pData(farms_informative2)$Treatment+pData(farms_informative2)$Animal, # BEFORE
        las= 2, cex.axis=0.7, main="BOLA-DQA2 = Bt.4751.1.S1_a_at (before)",ylim=c(6,14)) 
abline(h=means["Bt.4751.1.S1_a_at"])
animals.abline()
boxplot(tmp.eset["Bt.4751.1.S1_a_at",]~pData(farms2)$Treatment+pData(farms2)$Animal,                                          # AFTER
        las=2, cex.axis=0.7, main="BOLA-DQA2 = Bt.4751.1.S1_a_at (after)",ylim=c(6,14))
abline(h=means["Bt.4751.1.S1_a_at"])
animals.abline()
boxplot(exprs(farms_informative2)["Bt.25051.1.A1_at",]~pData(farms_informative2)$Treatment+pData(farms_informative2)$Animal,  # BEFORE
        las= 2, cex.axis=0.7, main="TTC14 = Bt.25051.1.A1_at (before)",ylim=c(6,14)) # boxplot of expression values for first of those probes
abline(h=means["Bt.25051.1.A1_at"])
animals.abline()
boxplot(tmp.eset["Bt.25051.1.A1_at",]~pData(farms2)$Treatment+pData(farms2)$Animal,                                          # AFTER
        las=2, cex.axis=0.7, main="TTC14 = Bt.25051.1.A1_at (after)",ylim=c(6,14))
abline(h=means["Bt.25051.1.A1_at"])
animals.abline()


# Let's look at a probe with no animal effect (esClas2 = Animal+Time)
which(esClas2[,"Int__Animal_724R.727R_vs_Treatment_BCG.CN"] != 0) # this line does the opposite,= looks for genes having an interaction effect involving animal ### Filter for genes DE in "Int__Animal_724R.727R_vs_Treatment_BCG.CN" column (+/- 1) in Animal+Timepoint

# I take all the contrasts involving "Animal" factor
tmp = esClas2[,grep("Animal", colnames(esClas2), ignore.case = FALSE, perl = FALSE, value = FALSE, # columns containing the "Animal" pattern
             fixed = FALSE, useBytes = FALSE, invert = FALSE)]

sum(rowSums(tmp) >= 1) # 3727 probes have at least an animal effect or an interaction effect involving animal
sum(rowSums(tmp) == 0) # 16437 probes have no significant animal effect

# We need to pick out an informative one (say a DE genes due to treatment)
# I will take one more useful column)
colSums(abs(esClas2)) # there are 94 DE genes between MAP and BOVIS (merging all time points)
# Let's use this one (721.MAP_vs_721.BOVIS )
tmp = esClas2[,grep("Animal", colnames(esClas2), ignore.case = FALSE, perl = FALSE, value = FALSE, # columns containing the "Animal" pattern
                    fixed = FALSE, useBytes = FALSE, invert = FALSE)]
tmp = cbind("721.MAP_vs_721.BOVIS"=esClas2[,"721.MAP_vs_721.BOVIS"], tmp)

rownames(tmp[(rowSums(abs(tmp)) == 1 & tmp[,"721.MAP_vs_721.BOVIS"] == 1),]) # 12 probes
# From the 94 genes which have DE between BOVIS and MAP all time points merged, 12 of them don't have any animal direct or interaction effect
#One of them is one described previously in this script: RAB7B / Bt.11413.1.S1_at
# Let's boxplot it (no animal effect means the values shouldb't change after correction)
par(mfrow=c(1,2))
boxplot(exprs(farms_informative2)["Bt.11413.1.S1_at",]~pData(farms_informative2)$Treatment+pData(farms_informative2)$Animal, # BEFORE
        las= 2, cex.axis=0.7, main="RAB7B = Bt.11413.1.S1_at (before)",ylim=c(6,14)) 
abline(h=means["Bt.11413.1.S1_at"])
animals.abline()
boxplot(tmp.eset["Bt.11413.1.S1_at",]~pData(farms2)$Treatment+pData(farms2)$Animal,                                          # AFTER
        las=2, cex.axis=0.7, main="RAB7B = Bt.11413.1.S1_at (after)",ylim=c(6,14))
abline(h=means["Bt.11413.1.S1_at"])
animals.abline()
#
#
#
#
#
# Top 9 BOLA probes = Top animal affected & annotated genes
# list of the 10 BOLA genes on top of the significance list of animal-influenced genes
tmp = list(
  probes=c("Bt.4751.1.S1_a_at","Bt.4751.2.S1_a_at","Bt.7220.1.S1_at","Bt.22867.1.S1_at","Bt.350.1.S1_x_at","Bt.4594.1.S1_at","Bt.350.1.S1_at",
    "Bt.29815.1.S1_x_at","Bt.22867.1.S1_x_at","Bt.29815.1.A1_at"),
  gene.symbols=c("BOLA-DQA2","BOLA-DQA2","BLA-DQB","BOLA-DQA1","BLA-DQB","BLA-DQB","BLA-DQB","BOLA","BOLA-DQA1","BOLA"))
par(mfrow=c(3,3))
# normalised uncorrected
for(i in 1:9){
  boxplot(assayData(farms2)$exprs[tmp$probes[i],]~pData(farms2)$Treatment+pData(farms2)$Animal, # BEFORE
          las= 2, cex.axis=0.7, main=paste(tmp$gene.symbols[i],"=",tmp$probes[i],"(uncorrected)",sep=" "),ylim=c(5,15))
  abline(h=means[tmp$probes[i]])
  animals.abline()
}

#corrected
for(i in 1:9){
  boxplot(assayData(farms2.animal.corrected)$exprs[tmp$probes[i],]~pData(farms2.animal.corrected)$Treatment+pData(farms2.animal.corrected)$Animal, # BEFORE
          las= 2, cex.axis=0.7, main=paste(tmp$gene.symbols[i],"=",tmp$probes[i],"(corrected)",sep=" "),ylim=c(5,15)) 
  abline(h=means[tmp$probes[i]])
  animals.abline()
}


##############
# Informative probes in the corrected dataset
##############

INIs <- INIcalls(farms2.animal.corrected) # Find informative genes from the farms normalised genes (probes?) using INIcalls
farms.informative.animal.corrected <- getI_Eset(INIs) # List of informative probesets, as another exprSet


##############
# Heatplot free of animal effect
##############

heatplot(dataset=farms.informative.animal.corrected, dualScale=FALSE, # faster for a subset of genes.
         scale="none", cexCol=0.5)
# Flipping the axis for better view
heatplot(dataset=t(assayData(farms.informative.animal.corrected)$exprs), dualScale=FALSE, # faster for a subset of genes.
         scale="none", cexCol=0.5, cexRow = 0.6)

# Increasing the PNG file width and heigth, to increase the axis labels to a readable size
png(filename="animal-corrected-heatplot-informative-huge.png", width=21, height=29.7, units="cm")
#heatplot(dataset=farms.informative.animal.corrected, dualScale=FALSE, # faster for a subset of genes.
#         scale="none", cexCol=0.5)
# For some reason, the heatplot above did not output anything after even a dozen minutes. 
# It could be that PNG files of A4 size take a long time to write, but I am not willing to wait hours to verify it.
dev.off()

#############
# PCA of corrected values
###########
# Hypothesis: do the samples cluster better by animal and timepoint after animal correction?

pca.corrected <- prcomp(data.frame(t(assayData(farms.informative.animal.corrected)$exprs)), scale=T) # clustering based on expression values
#summary(pca.corrected) # 25 to 5 % of the variance left after correction by each of the first 5 primary components
#summary(pca.corrected)$importance[, 1:6]  # The first 6 primary components explain 74% of the variance left after correction
#png(filename="PCA_normalised.png", width = 600, height = 600)
custom_col = c(rep("black",15),rep("red",15),rep("green",15),rep("blue",15)) # color by treatment (black, red, green, blue) for (CN, BCG, MAP, BOVIS), colors sorted blackRGB and treatments by virulence
custom_pch = rep(c(rep(3,5), rep(2,5),rep(7,5)),4) #3/2/7 shapes by time (+, A, *) for (2,6,24) The more lines in the symbol, the larger the time point
plot(pca.corrected$x, col=custom_col, pch=custom_pch,
     main="PCA - Corrected Probe Sets") 
# Lines below were copy-pasted from another PCA above in this script, adapt lines before running them
#legend(x="bottomright",col=legend_col, pch=legend_pch, inset = 0.01, cex=0.60,
#       legend=c("CN_0HR","CN_2HR","CN_6HR","CN_24HR","BCG_2HR","BCG_6HR","BCG_24HR",
#               "MAP_2HR", "MAP_6HR", "MAP_24HR", "BOV_2HR", "BOV_6HR","BOV_24HR"))
dev.off()


###########
# PCA of correction values (animal distance from the overall mean for each probe)
############
# Hypothesis: How many PC explain a fair % of (variance due to) the animal effect?

pca.animal.distance <- prcomp(data.frame(t(animals.distances), scale=T))
summary(pca.animal.distance) # Out of the 5 PC, only 4 are needed to explain 100% of the variance, 3 explain 84% of it
summary(pca.animal.distance)$importance
custom_col = 
custom_pch = 
plot(pca.animal.distance$x, #col=custom_col, pch=custom_pch,
     main="PCA - Animal distance") 

# Plotting the probes in animal space, which is a transformed matrix. Not sure how to interprete this transformation of the data

#pca.animal.distance2 <- prcomp(data.frame(animals.distances, scale=T))
#summary(pca.animal.distance2) # Out of the 5 PC, only 4 are needed to explain 100% of the variance, 3 explain 84% of it
#summary(pca.animal.distance2)$importance
#custom_col = 
#custom_pch = 
#plot(pca.animal.distance2$x, #col=custom_col, pch=custom_pch,
#       main="PCA - Animal distance") 


AnnotateProbeIDs = function(inlist)
{
  df = data.frame(ID=inlist)
  df$ID = inlist
  return(AnnotateTopTable(df))
}

tmp = names(which(abs(pca.animal.distance2$x[,3])>2))
AnnotateProbeIDs(tmp)

#############
# Which genes are driving the PC1, which explains 50% of the animal effect variance?
############
setwd("C:/Users/krue/Documents/Kevin-Logs/20130320_MDM/based-on-teles 2013")
tmp = names(which(abs(pca.animal.distance2$x[,1])>1))
write.table(x=AnnotateProbeIDs(tmp), file="351Probes-driving-PC1_sup1_animal-effect.txt", append=F, quote=F, sep="\t",
            eol="\n", dec=".", row.names=F, col.names=T)



###############
# Clean figures for Supervisor meeting (and potential publication)
###############
setwd("C:\\Users\\krue\\Documents\\Kevin-Logs\\20130320_MDM\\based-on-teles 2013\\For supervisor meeting") # Directory for the current project

################### PCA of uncorrected data
# Black, Red, Green, Blue code for CN, BCG, MAP, BOVIS
# Darker colors for larger time points
# Individually labelled animals
# Only informative probes used
opar = par()
par(mfrow=c(1,1)) # Make sure all the device is for just one plot
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE) # Add extra space to right of plot area; change clipping to figure
library(calibrate)

pca <- prcomp(data.frame(t(assayData(farms_informative2)$exprs)), scale=T) # clustering based on expression values
summary(pca)$importance[, 1:6]

data.frame(pca=rownames(pca$x), fi2=rownames(pData(farms_informative2))) # Overview of data rows to design correct coloring

colors = rep(NA, nrow(pData(farms_informative2)))
colors[which(pData(farms_informative2)$Treatment == "CN" & pData(farms_informative2)$TimePoint == "2HR")] = "grey80"
colors[which(pData(farms_informative2)$Treatment == "CN" & pData(farms_informative2)$TimePoint == "6HR")] = "grey40"
colors[which(pData(farms_informative2)$Treatment == "CN" & pData(farms_informative2)$TimePoint == "24HR")] = "black"
colfunc <- colorRampPalette(c("red3","white"))
colfunc(4)
colors[which(pData(farms_informative2)$Treatment == "BCG" & pData(farms_informative2)$TimePoint == "2HR")] = "#FFAAAA"
colors[which(pData(farms_informative2)$Treatment == "BCG" & pData(farms_informative2)$TimePoint == "6HR")] = "#FF5555"
colors[which(pData(farms_informative2)$Treatment == "BCG" & pData(farms_informative2)$TimePoint == "24HR")] = "#FF0000"
colors[which(pData(farms_informative2)$Treatment == "MAP" & pData(farms_informative2)$TimePoint == "2HR")] = "olivedrab1"
colors[which(pData(farms_informative2)$Treatment == "MAP" & pData(farms_informative2)$TimePoint == "6HR")] = "olivedrab3"
colors[which(pData(farms_informative2)$Treatment == "MAP" & pData(farms_informative2)$TimePoint == "24HR")] = "olivedrab4"
colfunc <- colorRampPalette(c("royalblue3","white"))
colfunc(4)
colors[which(pData(farms_informative2)$Treatment == "BOVIS" & pData(farms_informative2)$TimePoint == "2HR")] = "#BDC9EE"
colors[which(pData(farms_informative2)$Treatment == "BOVIS" & pData(farms_informative2)$TimePoint == "6HR")] = "#7B94DD"
colors[which(pData(farms_informative2)$Treatment == "BOVIS" & pData(farms_informative2)$TimePoint == "24HR")] = "#3A5FCD"

shapes = rep(NA, nrow(pData(farms_informative2)))
shapes[which(pData(farms_informative2)$TimePoint == "2HR")] = 17
shapes[which(pData(farms_informative2)$TimePoint == "6HR")] = 15
shapes[which(pData(farms_informative2)$TimePoint == "24HR")] = 18

size = rep(NA, nrow(pData(farms_informative2)))
size[which(pData(farms_informative2)$TimePoint == "2HR" | pData(farms_informative2)$TimePoint == "6HR")] = 1.5
size[which(pData(farms_informative2)$TimePoint == "24HR")] = 2

groups_legend = c("CN, 2H", "CN, 6H","CN, 24H", "BCG, 2H", "BCG, 6H","BCG, 24H", "MAP, 2H", "MAP, 6H","MAP, 24H", "BOVIS, 2H", "BOVIS, 6H","BOVIS, 24H")
colors_legend = c("grey80","grey40","black","#FFAAAA","#FF5555","#FF0000","olivedrab1","olivedrab3","olivedrab4","#BDC9EE","#7B94DD","#3A5FCD")
pch_legend = rep(c(17, 15, 18), 4)

plot(pca$x,  main="PCA - Informative Probe Sets", col=colors, pch=shapes, cex=size)  # PC1 vs PC2
textxy(pca$x[,1], pca$x[,2], pData(farms_informative2)$Animal,cx=.8)
legend("right", inset=c(-0.12,0), legend=groups_legend, pch=pch_legend, col=colors_legend, cex=1) # Add legend to top right, outside plot region
#plot(pca$x[,4:3],  main="PCA - Normalised Probe Sets", col=colors, pch=shapes, cex=size) # PC3 vs PC4
# For a fair comparison, I should compare PC1/PC2 before and after correction
# This shows that the treatment effect is the major effect, when the animal-biased probes are corrected (better than removed, which requires a threshold/criterium for removal)
plot(pca$x[,4:3],  main="PCA - Informative Probe Sets", col=colors, pch=shapes, cex=size)  # PC1 vs PC2
textxy(pca$x[,4], pca$x[,3], pData(farms_informative2)$Animal,cx=.8)
legend("right", inset=c(-0.12,0), legend=groups_legend, pch=pch_legend, col=colors_legend, cex=1) # Add legend to top right, outside plot region


################### PCA of corrected data

pca <- prcomp(data.frame(t(assayData(farms.informative.animal.corrected)$exprs)), scale=T) # clustering based on expression values
summary(pca)$importance[, 1:6]

plot(pca$x,  main="PCA - Informative corrected Probe Sets", col=colors, pch=shapes, cex=size)  # PC1 vs PC2
textxy(pca$x[,1], pca$x[,2], pData(farms_informative2)$Animal,cx=.8)
legend("right", inset=c(-0.12,0), legend=groups_legend, pch=pch_legend, col=colors_legend, cex=1) # Add legend to top right, outside plot region
#plot(pca$x[,4:3],  main="PCA - Normalised Probe Sets", col=colors, pch=shapes, cex=size) # PC3 vs PC4
# For a fair comparison, I should compare PC1/PC2 before and after correction
# This shows that the treatment effect is the major effect, when the animal-biased probes are corrected (better than removed, which requires a threshold/criterium for removal)
plot(pca$x[,4:3],  main="PCA - Normalised Probe Sets", col=colors, pch=shapes, cex=size)  # PC1 vs PC2
textxy(pca$x[,4], pca$x[,3], pData(farms_informative2)$Animal,cx=.8)
legend("right", inset=c(-0.12,0), legend=groups_legend, pch=pch_legend, col=colors_legend, cex=1) # Add legend to top right, outside plot region

par(opar) # revert to normal



################### Boxplots of expression
# Top BOLA, Last significant (TTC14), Not significant but DE between treatments (RAB7B)
par(mfrow=c(3,2))

uncorrected.boxplots = function(probe, gene)
{
  boxplot(exprs(farms_informative2)[probe,]~pData(farms_informative2)$Treatment+pData(farms_informative2)$Animal, # BEFORE
          las= 2, cex.axis=0.7, main=paste(gene,"=",probe,"(uncorrected)", sep=" "), ylim=c(6,14)) 
  abline(h=means[probe])
  animals.abline()
  animals.meanline(probe)
}

corrected.boxplots = function(probe, gene)
{
  boxplot(exprs(farms.informative.animal.corrected)[probe,]~pData(farms.informative.animal.corrected)$Treatment+pData(farms.informative.animal.corrected)$Animal, # BEFORE
          las= 2, cex.axis=0.7, main=paste(gene,"=",probe,"(corrected)", sep=" "), ylim=c(6,14)) 
  abline(h=means[probe])
  animals.abline()
  animals.meanline(probe)
}

animals.meanline = function(probe)
{
  i = 0
  for(animal in colnames(animals.means))
  {
    lines(x=c(1+4*i, 4*(i+1)), y=rep(animals.means[probe,animal],2), col="red", lty=4, lwd=1.75)
    i = i+1
  }
}

uncorrected.boxplots("Bt.4751.1.S1_a_at","BOLA-DQA2")
corrected.boxplots("Bt.4751.1.S1_a_at","BOLA-DQA2")
uncorrected.boxplots("Bt.25051.1.A1_at","TTC14")
corrected.boxplots("Bt.25051.1.A1_at","TTC14")
uncorrected.boxplots("Bt.11413.1.S1_at","RAB7B")
corrected.boxplots("Bt.11413.1.S1_at","RAB7B")

# Top 9 BOLAs before correction

uncorrected.boxplots5.14 = function(probe, gene)
{
  boxplot(exprs(farms_informative2)[probe,]~pData(farms_informative2)$Treatment+pData(farms_informative2)$Animal, # BEFORE
          las= 2, cex.axis=0.7, main=paste(gene,"=",probe,"(uncorrected)", sep=" "), ylim=c(5,14)) 
  abline(h=means[probe])
  animals.abline()
  animals.meanline(probe)
}

corrected.boxplots5.14 = function(probe, gene)
{
  boxplot(exprs(farms.informative.animal.corrected)[probe,]~pData(farms.informative.animal.corrected)$Treatment+pData(farms.informative.animal.corrected)$Animal, # BEFORE
          las= 2, cex.axis=0.7, main=paste(gene,"=",probe,"(corrected)", sep=" "), ylim=c(5,14)) 
  abline(h=means[probe])
  animals.abline()
  animals.meanline(probe)
}

par(mfrow=c(3,3))
uncorrected.boxplots5.14("Bt.4751.1.S1_a_at","BOLA-DQA2")
uncorrected.boxplots5.14("Bt.4751.2.S1_a_at","BOLA-DQA2")
uncorrected.boxplots5.14("Bt.7220.1.S1_at","BOLA-DQB")
uncorrected.boxplots5.14("Bt.22867.1.S1_at","BOLA-DQA1")
uncorrected.boxplots5.14("Bt.350.1.S1_x_at","BLA-DQB")
uncorrected.boxplots5.14("Bt.4594.1.S1_at","BLA-DQB")
uncorrected.boxplots5.14("Bt.350.1.S1_at","BLA-DQB")
uncorrected.boxplots5.14("Bt.29815.1.S1_x_at","BOLA")
uncorrected.boxplots5.14("Bt.22867.1.S1_x_at","BOLA-DQA1")


# Top 9 BOLAs after correction
corrected.boxplots5.14("Bt.4751.1.S1_a_at","BOLA-DQA2")
corrected.boxplots5.14("Bt.4751.2.S1_a_at","BOLA-DQA2")
corrected.boxplots5.14("Bt.7220.1.S1_at","BOLA-DQB")
corrected.boxplots5.14("Bt.22867.1.S1_at","BOLA-DQA1")
corrected.boxplots5.14("Bt.350.1.S1_x_at","BLA-DQB")
corrected.boxplots5.14("Bt.4594.1.S1_at","BLA-DQB")
corrected.boxplots5.14("Bt.350.1.S1_at","BLA-DQB")
corrected.boxplots5.14("Bt.29815.1.S1_x_at","BOLA")
corrected.boxplots5.14("Bt.22867.1.S1_x_at","BOLA-DQA1")


##############
# Differential expression of the corrected data
##############

# Fails when calculating a consensus_correlation for paired analysis.
# Say: nevermind
# when it worked at the 6th attempt, which was right after saving the session (unrelated i know, but weird anyway)

##################
# Expression boxplots of corrected data
##################
plot.corrected.expression = function(probe, gene)
{
  groups.corrected = groups[2:length(groups)]
  # Plots the normalised expression values for each condition for a given
  # probe
  # creates a temporary copy of the expression data, annotated by condition
  tmp_df = data.frame(t(assayData(farms.informative.animal.corrected)$exprs)) #
  tmp_df$group = NA #
  tmp_df$group[grep("CN_2HR",rownames(tmp_df))] = "CN_2HR" #
  tmp_df$group[grep("CN_6HR",rownames(tmp_df))] = "CN_6HR" #
  tmp_df$group[grep("CN_24HR",rownames(tmp_df))] = "CN_24HR" #
  tmp_df$group[grep("CN_25HR",rownames(tmp_df))] = "CN_24HR" #
  
  tmp_df$group[grep("BCG_2HR",rownames(tmp_df))] = "BCG_2HR" #
  tmp_df$group[grep("BCG_6HR",rownames(tmp_df))] = "BCG_6HR" #
  tmp_df$group[grep("BCG_24HR",rownames(tmp_df))] = "BCG_24HR" #
  tmp_df$group[grep("BCG_25HR",rownames(tmp_df))] = "BCG_24HR" #
  
  tmp_df$group[grep("MAP_2HR",rownames(tmp_df))] = "MAP_2HR" #
  tmp_df$group[grep("MAP_6HR",rownames(tmp_df))] = "MAP_6HR" #
  tmp_df$group[grep("MAP_24HR",rownames(tmp_df))] = "MAP_24HR" #
  tmp_df$group[grep("MAP_25HR",rownames(tmp_df))] = "MAP_24HR" #
  
  tmp_df$group[grep("BOVIS_2HR",rownames(tmp_df))] = "BOVIS_2HR" #
  tmp_df$group[grep("BOVIS_6HR",rownames(tmp_df))] = "BOVIS_6HR" #
  tmp_df$group[grep("BOVIS_24HR",rownames(tmp_df))] = "BOVIS_24HR" #
  tmp_df$group[grep("BOVIS_25HR",rownames(tmp_df))] = "BOVIS_24HR" #
  
  tmp_df$group = as.factor(tmp_df$group) #
  tmp_df$group = factor(tmp_df$group, groups.corrected) #
  boxplot(tmp_df[,probe]~group, data=tmp_df, cex.axis=0.78, las=3, main=paste(gene, "=", probe, "(corrected)",
          col=c(rep("green",3), rep("blue",3), rep("brown",3), rep("red",3)))
}

png(filename="expression-boxplot-corrected_BOLA-DQA2.png", width=700, height=700)
plot.corrected.expression("Bt.4751.1.S1_a_at", "BOLA-DQA2")
dev.off()
png(filename="expression-boxplot-corrected_TTC14.png", width=700, height=700)
plot.corrected.expression("Bt.25051.1.A1_at", "TTC14")
dev.off()
png(filename="expression-boxplot-corrected_RAB7B.png", width=700, height=700)
plot.corrected.expression("Bt.11413.1.S1_at", "RAB7B")
dev.off()

##############
# Differential expression of the corrected data
##############
setwd("C:/Users/krue/Documents/Kevin-Logs/20130320_MDM/R5")

# As stated above, the array correlation cannot be calculated (fails) for corrected data.

# To be consistent in this chimeric analysis, I will use the array_correlation value calculated for the uncorrected data
# Bear in mind that I already copied the se.exprs from the original data as well. Only the intensities have been corrected.

groups <- unique(targets$Group) # Group animals by (treatment,timepoint)
f <- factor(targets$Group, levels = groups) # let f be the tuple (treatment,timepoint) of each array

# Design the matrix for the 65 arrays, with the columns :
# "BOVIS_2HR","BOVIS_6HR","BOVIS_24HR","MAP_2HR","MAP_6HR","MAP_24HR"
design <- model.matrix(~0 + f)
colnames(design) <- groups
#design
rm(f)

# Design the matrix listing the contrasts to compute
cntrts <- c("BOVIS_2HR-MAP_2HR","BOVIS_6HR-MAP_6HR","BOVIS_24HR-MAP_24HR", "BOVIS_2HR-BCG_2HR","BOVIS_6HR-BCG_6HR","BOVIS_24HR-BCG_24HR", "BOVIS_2HR-CN_2HR","BOVIS_6HR-CN_6HR","BOVIS_24HR-CN_24HR", "MAP_2HR-BCG_2HR","MAP_6HR-BCG_6HR","MAP_24HR-BCG_24HR", "MAP_2HR-CN_2HR","MAP_6HR-CN_6HR","MAP_24HR-CN_24HR", "BCG_2HR-CN_2HR","BCG_6HR-CN_6HR","BCG_24HR-CN_24HR") # Prepares the contrasts names
contrasts <- makeContrasts(contrasts=cntrts, levels=design) # Make contrast matrix showing Bovis vs Map at each time point
#contrasts

# NOTE: peculiarity for paired comparison
#unique(targets$Animal)
animal <- targets$Animal # Show Animal column in targets file and call it animal
#animal

# The following block is required to obtain the correlation value between samples from the same animal
# this value will be required to block the analysis considering non-independence between samples from each animal

##Find out what is the correlation between all arrays from the same animal. Call it array_correlation.###
# The correlation is calculated from the unfiltered list of probe sets
array_correlation <- duplicateCorrelation(farms, design=design, ndups=1, block=animal)

##Show array_correlation###
#array_correlation
#array_correlation$consensus

############!!!!!!!!!!!####### And now we use farms.informative.animal.corrected again !!!
groups <- unique(pData(farms.informative.animal.corrected)$Group) # Group animals by (treatment,timepoint)
f <- factor(pData(farms.informative.animal.corrected)$Group, levels = groups) # let f be the tuple (treatment,timepoint) of each array
design <- model.matrix(~0 + f)
colnames(design) <- groups
rm(f)
cntrts <- c("BOVIS_2HR-MAP_2HR","BOVIS_6HR-MAP_6HR","BOVIS_24HR-MAP_24HR", "BOVIS_2HR-BCG_2HR","BOVIS_6HR-BCG_6HR","BOVIS_24HR-BCG_24HR", "BOVIS_2HR-CN_2HR","BOVIS_6HR-CN_6HR","BOVIS_24HR-CN_24HR", "MAP_2HR-BCG_2HR","MAP_6HR-BCG_6HR","MAP_24HR-BCG_24HR", "MAP_2HR-CN_2HR","MAP_6HR-CN_6HR","MAP_24HR-CN_24HR", "BCG_2HR-CN_2HR","BCG_6HR-CN_6HR","BCG_24HR-CN_24HR") # Prepares the contrasts names
contrasts <- makeContrasts(contrasts=cntrts, levels=design) # Make contrast matrix showing Bovis vs Map at each time point
animal <- pData(farms.informative.animal.corrected)$Animal # Show Animal column in targets file and call it animal

###Fit a linear model to the FARMS normalized data, call it farms_fit###
farms_fit = lmFit(farms.informative.animal.corrected, design = design, ndups=1, cor=array_correlation$consensus, block=animal)

###Apply contrasts to the linear model farms_fit###
farms_fit_contrast <- contrasts.fit(farms_fit, contrasts)

filter <- rownames(farms_fit_contrast) %in% rownames(exprs(farms_informative)) # Filter farms_fit_contrast using the informative genes found earlier (farms_informative)
filtered <- farms_fit_contrast[filter,] # long array of TRUE and FALSE whether to keep the row or not
# the above contains only the contrasted rows for the informative probesets 
# Below: Given a series of related parameter estimates and standard errors, compute moderated t-statistics, moderated F-statistic, and log-odds of differential expression by empirical Bayes shrinkage of the standard errors towards a common value.
eb_filtered <- eBayes(filtered)
#eb_filtered # Result of the differential expression statistical test, containing log odds of differential expression, t statistics, p.value, ...

# Write DEgenes to an outfile
write.table(eb_filtered, file = "DEgenes-eb_filtered-corrected.txt", quote = FALSE, sep = "\t", eol = "\n", na = "NA", row.names = TRUE, col.names = TRUE) 
# differential expression information of informative probes stored in a file


#####################
# log2FC 0 / Q 0.05 #
#####################
# Summary of which probe sets is DE in which contrast (useful for subsetting and Venn diagram)
# Note:  it will only conmsider one threshold: Q-value 0.05 and no threshold on fold-change
# probe name
# TRUE/FALSE for each contrast
# sum of TRUE
# gene.symbol
# ENTREZ+ENSEMBL IDs
DE.count.corrected = data.frame(row.names=1:iProbes) # initialises the data frame
DE.count.corrected$ID = topTable(eb_filtered, coef=1, number=iProbes, adjust.method="BH")$ID # gets the IDs of the informative probe sets
# Progressively merges the ID table with information about each contrast
for(cntr in cntrts){
  topT = topTable(eb_filtered, coef= cntr, number=iProbes, adjust.method="BH")
  topT$test = topT$adj.P.Val < 0.05
  names(topT)[names(topT)=="test"] <- cntr
  DE.count.corrected = merge(x=DE.count.corrected, y=topT[c("ID",cntr)], by="ID", sort=TRUE)
}
# clear memory
rm(topT)
# Summary stat: number of contrasts where each probe set was DE
DE.count.corrected$ContrastSum = rowSums(DE.count.corrected[names(DE.count.corrected)[grep(pattern="-", x=names(DE.count.corrected))]])
# Gene Annotation
iDE_count.corrected = AnnotateTopTable(DE.count.corrected) # 11952 rows, with duplicates from the 11842 informative probe sets
# write it to a file
write.table(x=iDE_count.corrected, file="Informative_probes_corrected_DE_contrasts.csv", append=F, quote=F, sep=",", eol="\n", dec=".", row.names=F)


################
# Cleanup
####################
rm(DE.df, Positive_DE_count, chimeric.contM, chimeric.design, chimeric.esClas, chimeric.eSet, chimeric.fit, chimeric.fitC, chimeric.pData, cn.FCs_df, common.df, contM, df,
   distances2, esetAll, fitC1.df, fitC2.df, pdataAll, res1, tmp.annotFC, tmp.common, tmp.cor, tmp.df, tmp.downs, tmp.eset, tmp.filtered, tmp.merged, tmp.subtract, tmp.ups,
   tmp_probes, INIs, esClas, esClas1, esClas2, farms3, filter, filtered, fit, fitC1, fitC2, pca, pca.animal.distance, pca.animal.distance2, pca.corrected, pdata, test, test.animal,
   test.measures, test.time, test.treatment, tmp, tmp.name, tmp.names, tmp.names.inter, tmp.names1, tmp.names2, tmp.var1plus)

################
# Ups and Downs
################

# For each contrast:
## get the full list of informative probes, sorted by P-value
## annotate the ordered list of probes with EnsEMBL and ENTREZ IDs
## add an entry in a summary table, of genes up and down (p-value or 0.05 and 0.01 should be enough)

# cntrts contains the list of contrast names. Use it to loop through the results

DE.df = data.frame()
for(cntr in cntrts){
  DE = topTable(eb_filtered, coef= cntr, number=iProbes, adjust.method="BH", sort.by="P")
  DE.df = UpsAndDowns(DE, cntr, DE.df)
  DE_annotated = AnnotateTopTable(DE)
  write.table(DE_annotated, file = paste(cntr,".txt", sep=""), quote = FALSE, sep = "\t", eol = "\n", na = "NA", row.names = FALSE, col.names = TRUE)
}
rm(DE, cntr, DE_annotated)
rownames(DE.df) = cntrts
colnames(DE.df) = c("up.FC0p05","down.FC0p05","DE.FC0p05","nDE.FC0p05",
                    "up.FC5p05","down.FC5p05","DE.FC5p05","nDE.FC5p05",
                    "up.FC0p01","down.FC0p01","DE.FC0p01","nDE.FC0p01")

write.table(x=DE.df, file="DEgenes_corrected-comparison-threshold.txt", quote=FALSE, sep="\t", eol = "\n", na = "NA", row.names = TRUE, col.names = TRUE)
rm(DE.df)


#############
# Custers of corregulated genes
############
setwd("C:/Users/krue/Documents/Kevin-Logs/20130320_MDM/Corregulated-corrected-genes")

correct.heatplot = heatplot(dataset=t(assayData(farms.informative.animal.corrected)$exprs), dualScale=FALSE, # faster for a subset of genes.
                            scale="none", cexCol=0.5, cexRow = 0.6, returnSampleTree=TRUE)

## Getting the members of a cluster and manuipulating the tree
## (from help file)
# str(correct.heatplot) # this draws  the full dendrogram in the standard out (console)
class(correct.heatplot)
plot(correct.heatplot)


abline(h=0.8, col="red")
## Cut the tree at the height=1.0
clusters.thirty = lapply(cut(correct.heatplot,h=0.8)$lower, labels)

# Trying to find the cluster containing IL1B (arbritrarily chosen to train myself to do this)
# doesn't seem able to find it anywhere 
which(clusters.thirty == "Bt.4856.1.S1_at") 
which("Bt.4856.1.S1_at" %in% clusters.thirty)
which(rownames(assayData(farms.informative.animal.corrected)$exprs) == "Bt.4856.1.S1_at") # the probe is in the corrected dataset!
for (l in clusters.thirty){
  print("Bt.4856.1.S1_at" %in% clusters.thirty)
}
# Maybe I use the list structure wrongly, I'll export the data in a TXT file and search the probe the "old way"
lapply(clusters.thirty, write, "30clusters-probes.txt", append=TRUE, ncolumns=10000)

# Alright, I found that Bt.4856.1.S1_at is in cluster 25
length(clusters.thirty[[25]]) # there are 677 probe members in that cluster
#Annotating the probes with gene.symbols
df = data.frame(ID=clusters.thirty[[25]])
df$ID = clusters.thirty[[25]]
df = AnnotateTopTable(df)


################
# Uncorrected boxplot of a non-BOLA gene: ATP2B1 (Bt.23296.1.S1_at)
##################
uncorrected.boxplots5.15(probe="Bt.23296.1.S1_at", gene="ATP2B1")

##################
# Measure of animal effect: sum of squared distance from each animal to overall mean, by probe set
##################
f.sum.squared = function(x){
  return(sum(x^2))
}
#head(apply(X=animals.distances, MARGIN=1, FUN="f.sum.squared")) # works
animal.effect.squared.sum = matrix(data=apply(X=animals.distances, MARGIN=1, FUN="f.sum.squared"),
                                   nrow=nrow(animals.distances), ncol=1,
                                   dimnames=list(rownames(animals.distances), "sq.sum.distances"))
# Let's check out BOLA-DQA2 = Bt.4751.1.S1_a_at
animal.effect.squared.sum["Bt.4751.1.S1_a_at",] # 41, compared to 0.033 for the first AFFX-BioB probes

# Let's rank and annotate this matrix
animal.effect.squared.sum = data.frame(animal.effect.squared.sum)
animal.effect.squared.sum$ID = rownames(animal.effect.squared.sum)
animal.effect.squared.sum = animal.effect.squared.sum[order(x=animal.effect.squared.sum$sq.sum.distances, decreasing=TRUE),]
animal.effect.squared.sum = AnnotateTopTable(animal.effect.squared.sum)

# Uncorrected boxplot of another non-BOLA gene top listed in this new list (suprisingly not in the PUMA calls)
uncorrected.boxplots5.15(probe="Bt.510.1.S1_at", gene="TAP")
uncorrected.boxplots5.15(probe="Bt.8124.1.S2_at",gene="COL1A2")
uncorrected.boxplots5.15(probe="Bt.2892.1.S1_at",gene="FABP7")

##################
# Trying to fit a GLM to COL1A2 (should have only an animal 724R positive effect, by eye)
##################
COL1A2.exprs = assayData(farms.informative.animal.corrected)$exprs["Bt.8124.1.S2_at",] # Getting the corrected expression values for COL1A2
animal_ = as.factor(pData(farms.informative.animal.corrected)$Animal) # Extracting group informations (x3)
timepoint_ = as.factor(pData(farms.informative.animal.corrected)$TimePoint)
treatment_ = as.factor(pData(farms.informative.animal.corrected)$Treatment)
COL1A2.glm.one = glm(COL1A2.exprs ~ 1 + animal_ + timepoint_ + treatment_, family=gaussian()) # Fitting an additive linear model, with gaussian distributed residuals
anova(COL1A2.glm.one)
summary(COL1A2.glm.one)
# OOOPS!! I run the analysis on the corrected intensities!!
# However, it called a statistically significant effect of the BOVIS treatment (p 0.05)
# Let's verify this by eye
corrected.boxplots5.14(probe="Bt.8124.1.S2_at", gene="COL1A2")
# The significant BOVIS effect is likely due to animal 724R corrected, where the BOVIS box(plot) does not intersect with the others

# Let's repeat the analysis as intended initially, to the uncorrected intensities
COL1A2.exprs = assayData(farms_informative2)$exprs["Bt.8124.1.S2_at",] # Getting the uncorrected expression values for COL1A2
animal_ = as.factor(pData(farms_informative2)$Animal) # Extracting group informations (x3)
timepoint_ = as.factor(pData(farms_informative2)$TimePoint)
treatment_ = as.factor(pData(farms_informative2)$Treatment)
COL1A2.glm.one.uncorrected = glm(COL1A2.exprs ~ 1 + animal_ + timepoint_ + treatment_, family=gaussian()) # Fitting an additive linear model, with gaussian distributed residuals
anova(COL1A2.glm.one.uncorrected)
summary(COL1A2.glm.one.uncorrected)
# YEAH !! THe linear model both captures the animal effect of animal 724R (+4.25891 log2 intensity) and the BOVIS effect seen after stripping out the animal effect.
# How biologically relevant this is, I am not sure though.

# Another thing which worries me
# is that the first level of each factor is not present in the summary table of the linear model
# I know this has to do with the degrees of freedom for each factor
# but I wonder how will we be able to identify an effect related to those levels of the factors?

# I need a probeset for which animal 706 is the only one standing apart from the others, say ... (?)
# I am lazy to chase after such a gene (which may or may not exist)
# Meanwhile, I have taken an example of BOLA-DQB 
# Where 706 and 724R cluster together, while 716,721, and 727R cluster together separately
DQB.exprs = assayData(farms_informative2)$exprs["Bt.7220.1.S1_at",] # Getting the uncorrected expression values for COL1A2
animal_ = as.factor(pData(farms_informative2)$Animal) # Extracting group informations (x3)
timepoint_ = as.factor(pData(farms_informative2)$TimePoint)
treatment_ = as.factor(pData(farms_informative2)$Treatment)
DQB.glm.one.uncorrected = glm(DQB.exprs ~ 1 + animal_ + timepoint_ + treatment_, family=gaussian()) # Fitting an additive linear model, with gaussian distributed residuals
anova(DQB.glm.one.uncorrected)
summary(DQB.glm.one.uncorrected)
# Here, 716, 721, and 727R are presented as having a significant animal effect
# Consequently, I expect that in a situations where 706 is the only outlier, in fact,
# the 4 other animals will be called as having a significant animal effect, relative to the 706 animal

# In short, animal 706, timepoint 2HR, and treatment BCG are used as references, representing the reduction of 1 degree of freedom.

# Now the best thing is to order the grouping factors, so that the references used are not automatically chosen as the first in alphabetical order
# but rather as meaningful references:
# - time: let's keep 2HR (0HR would have been the best if there were samples exposed to BCG/MAP/BOVIS and immediately measured)
# - treatment: CN
# - animal: well.. let's keep 704 as an arbritrary reference
animal_ = as.factor(pData(farms_informative2)$Animal) # Extracting group informations (x3)
timepoint_ = factor(pData(farms_informative2)$TimePoint, levels=c("2HR","6HR","24HR", ordered=TRUE))
treatment_ = factor(pData(farms_informative2)$Treatment,  levels=c("CN","BCG","MAP","BOVIS", ordered=TRUE))

# OK, let's try again BOLA-DQB
DQB.exprs = assayData(farms_informative2)$exprs["Bt.7220.1.S1_at",] # Getting the uncorrected expression values for COL1A2
DQB.glm.one.uncorrected = glm(DQB.exprs ~ 1 + animal_ + timepoint_ + treatment_, family=gaussian()) # Fitting an additive linear model, with gaussian distributed residuals
anova(DQB.glm.one.uncorrected)
summary(DQB.glm.one.uncorrected)
# YEAH!! 
# the animal effect is observed as animals 716, 721, and 727R being very different from animal 706 (when merging on treatment and time, those animals are very different from 706)
# an additional BCG and BOVIS effect is observed compared to CN (when merging on time and animals, BCG and BOVIS are sign. diff from CN)
# the time factor does not reach significance (when merging on treatment and anomal, there is no sign. difference between timepoints)

# In conclusion, this probeset is influenced by the animal and the treatment factor.
# Note: Because we don't have biological replicates for each factor combination, we have not enough data to assess any interaction effect.
uncorrected.boxplots5.15(probe="Bt.7220.1.S1_at",gene="BOLA-DQB")
plot.corrected.expression(probe="Bt.7220.1.S1_at",gene="BOLA-DQB")

###############
# Venn diagrams of corrected and uncorrected DE genes
##############
setwd("C:/Users/krue/Documents/Kevin-Logs/20130320_MDM/R5_corrected")
# We expect ALL the uncorrected DE genes to be included in the corrected list
# while the new list should include additional DE genes which treatment effect was overwhelmed by the animal effect
#colnames(DE.count)[3:ncol(DE.count) - 1] # no idea why I need to write 3 when I mean 2, but that is how it works
#sum(rownames(DE.count)==rownames(DE.count.corrected)) # rows of both datasets are identically ordered
par(mfrow=c(6,3))
for (contrast in colnames(DE.count)[3:ncol(DE.count) - 1]){
  vennDiagram(cbind(uncorrected=DE.count[,contrast], corrected=DE.count.corrected[,contrast]), main=contrast)
}
par = opar
#vennDiagram(cbind(uncorrected=DE.count[,"BOVIS_6HR-MAP_6HR"], corrected=DE.count.corrected[,"BOVIS_6HR-MAP_6HR"]), main="test")

##########
# Examples of gene DE after correction but not before
##########
DE.after.correction = function(contrast){
  return(DE.count[which(DE.count[,contrast] != DE.count.corrected[,contrast]),]$ID)
}

plot.DE.after.correction = function(contrast){
  opar = par()
  par(mfrow=c(2,2))
  print(contrast)
  print(sample(x=DE.after.correction(contrast), size=2, replace=FALSE))
  for(probeset in sample(x=DE.after.correction(contrast), size=2, replace=FALSE)){
    print(probeset)
    uncorrected.boxplots5.14(probe=probeset, gene=getSYMBOL(x=probeset, data="bovine.db"))
    corrected.boxplots5.14(probe=probeset, gene=getSYMBOL(x=probeset, data="bovine.db"))
  }
  par = opar
}

plot.DE.after.correction("BOVIS_6HR-MAP_6HR")

## RI ran the above twice, and obtained each time 2 randomly chosen probes 
# being DE between BOVIS and MAP at 6 hpi only after correction but not before

# The 4 genes involved: KLHL11, GBAS, ZNF48, CSF1R
opar = par()
par(mfrow=c(2,2))
plot_expression(probe="Bt.21907.1.S1_at", gene="ZNF48")
plot.corrected.expression(probe="Bt.21907.1.S1_at", gene="ZNF48")
plot_expression(probe="Bt.6151.1.S1_at", gene="CSF1R")
plot.corrected.expression(probe="Bt.6151.1.S1_at", gene="CSF1R")
par = opar

### Then I looked specifically for BOLA probe sets found DE after but not before correction
# Bt.22867.1.S1_x_at alias BOLA-DQA1 fits these criteria
opar = par()
par(mfrow=c(1,2))
uncorrected.boxplots5.14(probe="Bt.22867.1.S1_x_at", gene="BOLA-DQA1")
corrected.boxplots5.14(probe="Bt.22867.1.S1_x_at", gene="BOLA-DQA1")
par = opar

opar = par()
par(mfrow=c(1,2))
plot_expression(probe="Bt.22867.1.S1_x_at", gene="BOLA-DQA1")
plot.corrected.expression(probe="Bt.22867.1.S1_x_at", gene="BOLA-DQA1")
par = opar

opar = par()
par(mfrow=c(1,2))
# How fast does the animal effect between 706 and 721 weaken when going going down the ranking by p-value?
uncorrected.boxplots5.14(probe="Bt.29822.1.A1_at", gene="rank1")
uncorrected.boxplots5.14(probe="Bt.6432.1.S1_at", gene="rank400")
par = opar

#############
# IPA miRNA analysis of common DE genes at 2 and 6 hpi
############
# first define if each informative probe is "common at 2hrs" 
common.df = data.frame(rownames=DE.count$ID)
common.df$ID = DE.count$ID
common.df$common2 = apply(DE.count[,c("MAP_2HR-CN_2HR", "BCG_2HR-CN_2HR", "BOVIS_2HR-CN_2HR")], MARGIN=1, FUN=all)
common.df$common6 = apply(DE.count[,c("MAP_6HR-CN_6HR", "BCG_6HR-CN_6HR", "BOVIS_6HR-CN_6HR")], MARGIN=1, FUN=all)
common.df$common24 = apply(DE.count[,c("MAP_24HR-CN_24HR", "BCG_24HR-CN_24HR", "BOVIS_24HR-CN_24HR")], MARGIN=1, FUN=all)
# 2 hours
setwd("C:/Users/krue/Documents/Kevin-Logs/20130320_MDM/IPA-Overlap/")
tmp.common = AnnotateTopTable(common.df)
head(tmp.common[tmp.common$common2,])
unique(tmp.common[tmp.common$common2,]$ID) # 415 (unique probe sets)
write.table(unique(tmp.common[tmp.common$common2,]$ID), file="probesets_common_2hpi.txt", quote=F, sep="\t", eol="\n", dec=".", row.names=F, col.names=F)
#common.df[common.df$common2,]$ID # 415
#write.table(common.df[common.df$common2,]$ID, file="probeIDs_common_2hpi.txt", quote=F, sep="\t", eol="\n", dec=".", row.names=F, col.names=F)
rm(tmp.common)
# 6 hours
tmp.common = AnnotateTopTable(common.df)
head(tmp.common[tmp.common$common6,])
unique(tmp.common[tmp.common$common6,]$ID) # 415 (unique probe sets)
write.table(unique(tmp.common[tmp.common$common6,]$ID), file="probesets_common_6hpi.txt", quote=F, sep="\t", eol="\n", dec=".", row.names=F, col.names=F)
#common.df[common.df$common2,]$ID # 415
#write.table(common.df[common.df$common2,]$ID, file="probeIDs_common_2hpi.txt", quote=F, sep="\t", eol="\n", dec=".", row.names=F, col.names=F)
rm(tmp.common, common.df)

rm(common.df)

# This does not work. IPA can work with a list of miRNA to find enriched targets
# but not the other way around. I have a gene list.


###################
# A few more animal boxplots
###################
# Papers about "inter-animal variability" mention cytochromes (CYP...) as variable between animals
# CYP39A1, CYP27A1, CYP4V2 are indeed DE between animal 706 and 721 according to my first PUMA analyses above
uncorrected.boxplots4.15(probe="Bt.14369.1.A1_at", "CYP39A1") # very clear differences between animals
uncorrected.boxplots4.15(probe="Bt.16001.1.S1_at", "CYP27A1") # weaker animal effect here an below
uncorrected.boxplots4.15(probe="Bt.18730.1.A1_at", "DHFR")
uncorrected.boxplots4.15(probe="Bt.13224.1.A1_at", "CYP4V2")
uncorrected.boxplots4.15(probe="Bt.16514.1.S1_at", "CYP4V2")
uncorrected.boxplots4.15(probe="Bt.9714.1.S1_at", "ACTG2")


###################
# Comparison MDM vs RNA-seq for top 100 probesets with animal effect
#####################
setwd("C:/Users/krue/Documents/Kevin-Logs/20130320_MDM/microarray-RNAseq-comparison")

# First, calculate a metric of the animal effect
# = sum of squared distance of all animals' mean to the overall mean for each probeset

f.sum.squared = function(x){
  return(sum(x^2))
}
#head(apply(X=animals.distances, MARGIN=1, FUN="f.sum.squared")) # works
animal.effect.squared.sum = matrix(data=apply(X=animals.distances, MARGIN=1, FUN="f.sum.squared"),
                                   nrow=nrow(animals.distances), ncol=1,
                                   dimnames=list(rownames(animals.distances), "sq.sum.distances"))
# Let's check out BOLA-DQA2 = Bt.4751.1.S1_a_at
animal.effect.squared.sum["Bt.4751.1.S1_a_at",] # 41, compared to 0.033 for the first AFFX-BioB probes

# Let's rank and annotate this matrix
animal.effect.squared.sum = data.frame(animal.effect.squared.sum)
animal.effect.squared.sum$ID = rownames(animal.effect.squared.sum)
animal.effect.squared.sum = animal.effect.squared.sum[order(x=animal.effect.squared.sum$sq.sum.distances, decreasing=TRUE),]
animal.effect.squared.sum = AnnotateProbesTable(df=animal.effect.squared.sum, annotPkg="bovine.db",
                                                probeCol="ID", gene.symbol=F, ENSEMBL=T, ENTREZ=F)

# Merge this table with the expression intensities (without animal correction)
subset1 = colnames(assayData(farms2)$exprs)[grep(pattern="BOVIS_2.HR", x=colnames(assayData(farms2)$exprs))]
subset2 = colnames(assayData(farms2)$exprs)[grep(pattern="CN_2.HR", x=colnames(assayData(farms2)$exprs))]
farms2.animalDistance.exprs = merge(x=animal.effect.squared.sum, y=assayData(farms2)$exprs[,c(subset1, subset2)], by.x="ID", by.y="row.names")
farms2.animalDistance.exprs = farms2.animalDistance.exprs[order(farms2.animalDistance.exprs$sq.sum.distances, decreasing=TRUE),]

# Takes only the informative probesets
farms2.animalDistance.exprs = farms2.animalDistance.exprs[which(farms2.animalDistance.exprs$ID %in% rownames(assayData(farms_informative2)$exprs)),]


# Take the first 100 rows of this dataset (ranked by inter-animal average distance)
save(farms2.animalDistance.exprs, file="farms2.animalDistance.exprs.RData")