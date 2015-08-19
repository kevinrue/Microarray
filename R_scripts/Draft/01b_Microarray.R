


# Libraries ---------------------------------------------------------------


library(affy)
library(arrayQualityMetrics)
library(farms)
library(ggplot2)
library(RColorBrewer)
library(reshape2)
library(tools)


# General set up ----------------------------------------------------------


## directory with CEL files
CEL.dir = 'C:/Users/krue/Documents/Kevin-Logs/CELfiles/Project/ProjectData/CEL Files detailed'


# Folder where the results of the analysis will be saved
rootDir = 'C:/Users/krue/Documents/Kevin-Logs/20130320_MDM/R6_complete_analysis'

setwd(rootDir)


# Import and format the phenotypic information ----------------------------


# lists all the files in the directory
phenodata = data.frame(Filename = dir(CEL.dir), stringsAsFactors=F)
head(phenodata)

# Uses the filename , withut ".CEL", as the ID
phenodata$ID <- unlist(strsplit(phenodata$Filename, split = '.CEL', fixed = TRUE))
head(phenodata)

# Creates a column for the Animal_ID, the Infection, and the timpoint
temp.list = unlist(strsplit(phenodata$ID, split = "_", fixed = FALSE))
phenodata$Animal <- factor(temp.list[seq(1, length(temp.list), 3)])
phenodata$Infection <- factor(
  temp.list[seq(2, length(temp.list), 3)],
  levels = c('CN', 'BCG', 'MAP', 'BOVIS'),
  labels = c('Control', 'M. bovis BCG', 'MAP', 'M. bovis')
  )
# Simultaneously replace the 25HR Timepoint value by 24HR
phenodata$Timepoint <- factor(
  gsub(
    pattern = '25HR',
    replacement = '24HR',
    x = temp.list[seq(3, length(temp.list), 3)]),
  levels = c('0HR', '2HR', '6HR', '24HR'),
  labels = c('0 hpi', '2 hpi', '6 hpi', '24 hpi'),
  ordered = TRUE
  )
rm(temp.list)
head(phenodata)

# Check that the new columnns are factors
phenodata$Animal
phenodata$Infection
phenodata$Timepoint

# Reorder the samples by time point, then infection, then animal
phenodata = phenodata[order(
  phenodata$Infection, phenodata$Timepoint, phenodata$Animal),
  ]

# Create a column for the experimental group
phenodata$Group = factor(
  paste(phenodata$Infection, phenodata$Timepoint, sep = ' '),
  levels = unique(paste(phenodata$Infection, phenodata$Timepoint, sep = ' '))
  )

# Add a column to contain the timepoint as a numeric value
phenodata$Hours = as.numeric(
  gsub(pattern = ' hpi', replacement = '', x = phenodata$Timepoint)
  )

# Export the phenodata to a text file
write.table(
  phenodata, file="phenodata.txt", append=FALSE, quote=FALSE, sep="\t",
  eol="\n", row.names=FALSE, col.names=TRUE
  )


# Import microarray data --------------------------------------------------


# Read all CEL files into an affybatch object
rawdata <- ReadAffy(
  filenames=phenodata$Filename,
  celfile.path=CEL.dir,
  phenoData=AnnotatedDataFrame(phenodata),
  sampleNames=phenodata$ID)
# Examine the affybatch
rawdata

# Create a folder to save the quality control files
dir.create(file.path(rootDir, 'QC'))

# Quality control
arrayQualityMetrics(
  expressionset = rawdata,
  outdir = file.path(rootDir, 'QC'),
  do.logtransform = TRUE
  )

# 700_BCG_2HR was found an outlier in 5 out of 6 tests
# 713_MAP_6HR has much higher RNA degradation and failed 2 out of 6 tests

# Discard all samples from animals 700 and 713 from the phenodata
phenodata.65 = phenodata[!phenodata$Animal %in% c('700', '713'),]

# Export the new phenodata to a text file
write.table(
  phenodata.65, file="phenodata.65.txt", append=FALSE, quote=FALSE, sep="\t",
  eol="\n", row.names=FALSE, col.names=TRUE
)

# Cleanup
rm(phenodata, rawdata)



# Import the good quality microarray data ---------------------------------


# Read all CEL files into an affybatch object
rawdata.65 <- ReadAffy(
  filenames=phenodata.65$Filename,
  celfile.path=CEL.dir,
  phenoData=AnnotatedDataFrame(phenodata.65),
  sampleNames=phenodata.65$ID)
# Examine the affybatch
rawdata.65

# Create a folder to save the quality control files
dir.create(file.path(rootDir, 'QC.65'))

# Quality control
arrayQualityMetrics(
  expressionset = rawdata.65,
  outdir = file.path(rootDir, 'QC.65'),
  do.logtransform = TRUE
)


# Normalisation -----------------------------------------------------------


# Normalize the object rawdata using the FARMS normalization method
farms <- qFarms(rawdata.65)
# Cleanup the object not needed anymore
rm(rawdata.65)

# Boxplot of normalised data
tmp.matrix = exprs(farms)
tmp.df = melt(tmp.matrix, value.name = 'Intensity', varnames = c('Probeset', 'Sample'))
pdf(file = file.path(rootDir, 'Normalised.boxplot.pdf'), width = 14)
ggplot(
  data = tmp.df,
  aes(
    x = Sample,
    y = Intensity)
  ) +
  geom_boxplot(
    fill = brewer.pal(3, 'Pastel1')[1]
    ) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1)
    ) +
  ylab(expression('Log'[2]*' intensity'))
dev.off()

