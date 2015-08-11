
#############
# Libraries #
#############

library(affy) # loads the package
library(arrayQualityMetrics) #
library(affyPLM) #
library(simpleaffy) #
library(farms) #
library(limma) #
library(bovine.db) #
library(annotate) #
library(puma) #
library(made4) #
library(ggplot2) #
library(PerformanceAnalytics)

################### It only includes 5 animals passing QC for all samples.
## Data import 4 ##
###################
setwd("C:\\Users\\krue\\Documents\\Kevin-Logs\\20130320_MDM\\R4_QC-filtered") # Directory for the current project

# read annotated file into targets information
targets = read.AnnotatedDataFrame ("phenodata.txt", header = TRUE, row.names = 2, as.is = TRUE)

# get sample names from targets
SampleNames <- sampleNames(targets)

setwd("C:\\Users\\krue\\Documents\\Kevin-Logs\\CELfiles\\Project\\ProjectData\\CEL Files detailed")

# Read all CEL files into an affybatch object called rawdata
rawdata <- ReadAffy(filenames=pData(targets)$FileName, celfile.path=".", phenoData=targets, sampleNames=SampleNames)

# Back to the working directory
setwd("C:\\Users\\krue\\Documents\\Kevin-Logs\\20130320_MDM\\R4_QC-filtered")

# Saves the rawdata variable if we need it later
save(rawdata, file="rawdata.RData")

# clears the workspace
rm(SampleNames, targets, rawdata)

#################
# Normalisation #
#################
# Normalize the object rawdata using the FARMS normalization method. Call the normalized object farms 
setwd("C:\\Users\\krue\\Documents\\Kevin-Logs\\20130320_MDM\\R4_QC-filtered")
load("rawdata.RData")
farms = qFarms(rawdata)
rm(rawdata)

# saves data and clears workspace
save(farms, file="farms.RData")
rm(farms)

#####################
# Informative probe #
#####################
setwd("C:\\Users\\krue\\Documents\\Kevin-Logs\\20130320_MDM\\R4") # Directory for the current project

# Informative probes
load("farms.RData")
INIs <- INIcalls(farms) # Find informative genes from the farms normalised genes (probes?) using INIcalls
rm(farms)

farms_informative <- getI_Eset(INIs) # List of informative probesets, as another exprSet

# saves data and clears workspace
save(farms_informative, file = "farms_informative.RData")
rm(INIs, farms_informative)
