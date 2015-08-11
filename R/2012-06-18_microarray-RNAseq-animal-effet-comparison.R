# This file is meant for the comparison of RNAseq and microarray measures of gene expression
# for genes with an animal effect

# Remember to load the TranscriptomicsFunctions.R script

# Rule: Load data when needed, and remove after used

setwd("C:/Users/krue/Documents/Kevin-Logs/20130320_MDM/microarray-RNAseq-comparison/")

# Load the normalised 11842 informative microarray probesets
load(file="C:/Users/krue/Documents/Kevin-Logs/20130320_MDM/R4_QC-filtered/farms_informative.RData")

# Subset microarray samples: informative, CN+BOVIS, 24HR (11842 probesets)
farms.inf.BovCn.24 = farms_informative[,
                  intersect(
                    x=which(pData(farms_informative)$TimePoint == "24HR"),
                    y=which(pData(farms_informative)$Treatment == "CN" | pData(farms_informative)$Treatment == "BOVIS"))
                  ]
# Keep the description of samples (valuable for boxplotting)
pdata.BovCn.24 = pData(farms.inf.BovCn.24)
save(pdata.BovCn.24, file="pdata.BovCn.24.RData")
rm(farms_informative)

# Microarray expression data: informative, CN+BOVIS, 24HR (11842 probesets)
microarry.exprs = data.frame(assayData(farms.inf.BovCn.24)$exprs)
colnames(microarry.exprs) = colnames(assayData(farms.inf.BovCn.24)$exprs)
microarry.exprs$ID = rownames(assayData(farms.inf.BovCn.24)$exprs)
save(microarry.exprs, file="microarry.exprs")
rm(farms.inf.BovCn.24)

# Annotated Microarray expression data: informative, CN+BOVIS, 24HR (11952 rows) (n:n relationship probeset : multiple ensembl_id)
microarray.ensembl = AnnotateProbesTable(df=microarry.exprs, annotPkg="bovine.db", probeCol="ID", gene.symbol=F, ENSEMBL=T, ENTREZ=F)
save(microarray.ensembl, file="microarray.ensembl.RData")
rm(microarry.exprs)

# There are 3132 probesets without ensembl annotations
sum(is.na(microarray.ensembl$ensembl_id))
# There are indeed 11842 unique probesets_ids in the dataset
length(unique(microarray.ensembl$ID))
# There are 6879 unique ensembl_ids in the dataset
length(unique(microarray.ensembl$ensembl_id))
# There are 0 row without a probeset_id
sum(is.na(microarray.ensembl$ID))

# Annotated Microarray expression data: informative, CN+BOVIS, 24HR, with ensembl annotation (8820 rows) (n:n relationship probeset : multiple ensembl_id)
microarray.ensembl.noNA = microarray.ensembl[!is.na(microarray.ensembl$ensembl_id),]
save(microarray.ensembl.noNA, file="microarray.ensembl.noNA.RData")
rm(microarray.ensembl, microarray.ensembl.noNA)

# All 8820 rows have a unique probeset/ensembl_id pair
nrow(unique(microarray.ensembl.noNA[,c("ID","ensembl_id")]))

# Subset of ensembl IDs matching the following conditions:
## microarray: informative, ensembl annotation
## RNAseq:     non-zero count in all samples
load("C:/Users/krue/Documents/Kevin-Logs/20130320_MDM/microarray-RNAseq-comparison/rpkm.df.noZero.RData")
load("C:/Users/krue/Documents/Kevin-Logs/20130320_MDM/microarray-RNAseq-comparison/microarray.ensembl.noNA.RData")
inter = intersect(rownames(x=rpkm.df.noZero), y=microarray.ensembl.noNA$ensembl_id)

microarray.rnaseq.intersection = microarray.ensembl.noNA[microarray.ensembl.noNA$ensembl_id %in% inter,]
save(microarray.rnaseq.intersection, file="microarray.rnaseq.intersection.RData")
rm(microarray.ensembl.noNA, microarray.rnaseq.intersection)

rpkm.microarray.intersection = rpkm.df.noZero[rownames(rpkm.df.noZero) %in% inter,]
rpkm.microarray.intersection = log2(rpkm.microarray.intersection)
save(rpkm.microarray.intersection, file="rpkm.microarray.intersection")
rm(rpkm.df.noZero, rpkm.microarray.intersection)
rm(inter)


# average distance for each animal (all treatments, all timepoints) from the meanof all 5 animals for each probe/ensembl_id
# Ranked by decreasing sum of squared distance from the mean (i.e. decreasing animal effect)
# Filtered for ensembl_id in the above intersection
#
# A more accurate metric would be to divide this sum by the average deviation within each animal (variance between / variance within)
sq.distance = farms2.animalDistance.exprs[,1:3] # See farms2.animalDistance.exprs in the master script 2012-02-20_MDM-all.R
sq.distance = sq.distance[sq.distance$ensembl_id %in% inter,]
save(sq.distance, file="sq.distance.RData")
rm(farms2.animalDistance.exprs, sq.distance)
load("sq.distance.RData")

# Boxplot function a given probe in the microarray dataset
boxplot_microarray.animal = function(ensembl, probeset)
{
  load("C:/Users/krue/Documents/Kevin-Logs/20130320_MDM/microarray-RNAseq-comparison/pdata.BovCn.24.RData")
  load("C:/Users/krue/Documents/Kevin-Logs/20130320_MDM/microarray-RNAseq-comparison/microarray.rnaseq.intersection.RData")
  # Renames for convenience
  df = microarray.rnaseq.intersection
  annot = pdata.BovCn.24
  # library loading
  library(annotate, quietly=TRUE)
  library(bovine.db, quietly=TRUE)
  # Picks the data corresponding to the ensembl
  data = data.frame(probe=t(df[df$ID == probeset & df$ensembl_id == ensembl,2:11]))
  colnames(data) = "probe"
  data$group = annot$Animal
  # Gene Symbol for clarity
  gene.symbol=getSYMBOL(probeset, "bovine.db")
  # Boxplot
  boxplot(probe~group, data=data, ylim=c(4, 16),
          main = paste("Microarray:\n",probeset," = ", gene.symbol,"\n",ensembl),
          ylab = "log2(Intensity)")
}


#boxplot_microarray.animal(ensembl="ENSBTAG00000012341", probeset="Bt.1.1.S1_at")


# Boxplot function a given probe in the RNAsq dataset
boxplot_rnaseq.animal = function(ensembl)
{
  load("C:/Users/krue/Documents/Kevin-Logs/20130320_MDM/microarray-RNAseq-comparison/rpkm.microarray.intersection")
  df = rpkm.microarray.intersection
  annot = c("706","706","716","716","721","721","724R","724R","727R","727R")
  # Picks the data 
  data = data.frame(probe=t(df[rownames(df) == ensembl,]))
  colnames(data) = "probe"
  data$group = annot
  # Boxplot
  boxplot(probe~group, data=data, ylim=c(-8, 14),
          main = paste("RNAseq: ", ensembl),
          ylab = "log2(RPKM)")
}

#boxplot_rnaseq.animal(ensembl="ENSBTAG00000019486")


boxplot_all = function(ensembl)
{
  # First checks how many probes exist for this ensembl ID
  load("C:/Users/krue/Documents/Kevin-Logs/20130320_MDM/microarray-RNAseq-comparison/microarray.rnaseq.intersection.RData")
  probesets = microarray.rnaseq.intersection[microarray.rnaseq.intersection$ensembl_id == ensembl,]$ID
  n = length(probesets)
  if (n == 0 ){
    stop("No data for ensembl ID: ", ensembl)
  }
  # Prepares the graphical window
  columns = ceiling(sqrt(n+1))
  rows = ceiling((n+1)/columns)
  par(mfrow=c(rows,columns))
  # Boxplots the RNAseq data
  boxplot_rnaseq.animal(ensembl=ensembl)
  # Boxplots the microarray data
  for (probeset in probesets){
    boxplot_microarray.animal(ensembl=ensembl, probeset=probeset)
  }
}


# There we go, the ranked list of annotated probesets by decreasing mean animal distance
# Although, there is a lot of redundancy at the top, the same BOLA genes are measured by multiple probesets
# I will take the top 10 unique
boxplot_all(ensembl=unique(sq.distance$ensembl_id)[1]) # 9656 /  BOLA-DQA1/BOLA-DQA2
## Interesting because animal 721 is also low in RNAseq, and always low in microarray ##
## While all other animals have at least one probeset high ##
boxplot_all(ensembl=unique(sq.distance$ensembl_id)[2]) # 19588 / BLA-DQB # RNAseq always high, while microarray variable
## Probeset 3500.1.S1_s has a profile similar to RNAseq ##
boxplot_all(ensembl=unique(sq.distance$ensembl_id)[3]) # 2069 / BOLA # RNAseq always high, while microarray variable
## Probeset 3500.1.S1_s has a profile similar to RNAseq ##
boxplot_all(ensembl=unique(sq.distance$ensembl_id)[4]) # 9552 / ATP2B1 # RNAseq always high, while microarray variable
boxplot_all(ensembl=unique(sq.distance$ensembl_id)[5]) # 48171 / TAP
## Interesting because of similar profile for all animals ##
boxplot_all(ensembl=unique(sq.distance$ensembl_id)[6]) # 2087 / IL36A
## Interesting because of similar profile for all animals ##
boxplot_all(ensembl=unique(sq.distance$ensembl_id)[7]) # 14402 / GIMAP8 # RNAseq always high, while microarray variable
## Although profiles of RNAseq and one probeset are fairly similar
boxplot_all(ensembl=unique(sq.distance$ensembl_id)[8]) # 9210 / ZBTB44 # RNAseq always high, while microarray consistently different between animals
boxplot_all(ensembl=unique(sq.distance$ensembl_id)[9]) # 11985 / LOC509034
## Interesting because of similar profile for all animals ##
boxplot_all(ensembl=unique(sq.distance$ensembl_id)[10]) # 4950 / BRB
## Interesting because of similar profile for all animals ##

# Least animal affected genes
boxplot_all(ensembl=sq.distance[nrow(sq.distance),]$ensembl_id) # 9051 / SRP19 # Constant across animals
boxplot_all(ensembl=sq.distance[nrow(sq.distance)-1,]$ensembl_id) # 3443 / USP46 # Constant across animals
boxplot_all(ensembl=sq.distance[nrow(sq.distance)-2,]$ensembl_id) # 12199 / IMP3 # Constant across animals
boxplot_all(ensembl=sq.distance[nrow(sq.distance)-3,]$ensembl_id) # 12199 / TPD52L2 # Constant across animals
boxplot_all(ensembl=sq.distance[nrow(sq.distance)-4,]$ensembl_id) # 9742 / PDCD6IP # Constant across animals


