#
# This script is meant to analyse the correlation of RNAseq and microarray
# for probesets which are filtered on their alignment to the genome
# The filter is based on the percentage of probes from a same probeset
# which aligned to the most represented chromosome for that probeset.
#

# The probesets + mode(chromosome) + % aligned to chromosome
load(file="modes_percentages_total.RData")

# log2(RPKM) for genes with both RNAseq and microarray informative data
load(file="rpkm.microarray.intersection.RData")

# log2(intensity) for genes with both RNAseq and microarray informative data
load(file="microarray.rnaseq.intersection.RData")

############### Fist gene on array with 100% probes mapping to same chromosome

# Probesets 100% aligned to same chromosome
p1 = modes.percentages[modes.percentages$perc == 1,]$probeset
# Filter for probesets interesting with RNAseq data
p1.inter = intersect(x=p1, y=microarray.rnaseq.intersection$ID)
# An example
p1.inter[1] # "Bt.1.1.S1_at"
# The related ensembl_id
microarray.rnaseq.intersection[microarray.rnaseq.intersection$ID == "Bt.1.1.S1_at","ensembl_id"] # "ENSBTAG00000012341"

# Plot the microarray and RNAseq data for that probeset
order.micro = c("706_CN_24HR","716_CN_24HR","721_CN_24HR","724R_CN_25HR","727R_CN_25HR", 
              "706_BOVIS_24HR","716_BOVIS_24HR","721_BOVIS_24HR","724R_BOVIS_25HR","727R_BOVIS_25HR")
order.rna = c("706cn25","716cn25","721cn25","724Rcn25","727Rcn25","706bovis25","716bovis25","721bovis25","724Rbovis25","727Rbovis25")

plot(x = as.numeric(rpkm.microarray.intersection[rownames(rpkm.microarray.intersection) == "ENSBTAG00000012341", order.rna]),
     y = as.numeric(microarray.rnaseq.intersection[microarray.rnaseq.intersection == "Bt.1.1.S1_at",order.micro]),
     xlab = "RNAseq",
     ylab = "Microarray",
     main = paste("Bt.1.1.S1_at = ENSBTAG00000012341\n",
                  round(modes.percentages[modes.percentages$probeset == "Bt.1.1.S1_at",]$perc * 100, digits=2), "% of ",
                  modes.percentages[modes.percentages$probeset == "Bt.1.1.S1_at",]$total, " probes aligned", sep=""))

cor(x = as.numeric(rpkm.microarray.intersection[rownames(rpkm.microarray.intersection) == "ENSBTAG00000012341", order.rna]),
    y = as.numeric(microarray.rnaseq.intersection[microarray.rnaseq.intersection == "Bt.1.1.S1_at",order.micro]))

# Hmmm... well, poor negative correlation. This gene has little absolute variation anyway. 

############## BOLA DQA2

# We know that BOLA-DQA2 correlates visually very well between probeset Bt.4751.2.S1_a_at and RNAseq for ENSBTAG00000009656
microarray.rnaseq.intersection[microarray.rnaseq.intersection$ID == "Bt.4751.2.S1_a_at","ensembl_id"] # "ENSBTAG00000009656"
# animal 721 is an outlier in both technologies
plot(x = as.numeric(rpkm.microarray.intersection[rownames(rpkm.microarray.intersection) == "ENSBTAG00000009656", order.rna]),
     y = as.numeric(microarray.rnaseq.intersection[microarray.rnaseq.intersection == "Bt.4751.2.S1_a_at",order.micro]),
     xlab = "RNAseq",
     ylab = "Microarray",
     main = paste("Bt.4751.2.S1_a_at = ENSBTAG00000009656\n",
                  round(modes.percentages[modes.percentages$probeset == "Bt.4751.2.S1_a_at",]$perc * 100, digits=2), "% of ",
                  modes.percentages[modes.percentages$probeset == "Bt.4751.2.S1_a_at",]$total, " probes aligned", sep=""))
# cor 0.98
cor(x = as.numeric(rpkm.microarray.intersection[rownames(rpkm.microarray.intersection) == "ENSBTAG00000009656", order.rna]),
    y = as.numeric(microarray.rnaseq.intersection[microarray.rnaseq.intersection == "Bt.4751.2.S1_a_at",order.micro]))


########## General analysis
cor.RNA.micro = function(x)
{
  # Clear processing
  probeset = x["ID"]
  ensembl = x["ensembl_id"]
  # Orders the samples for paired correlation
  order.micro = c("706_CN_24HR","716_CN_24HR","721_CN_24HR","724R_CN_25HR","727R_CN_25HR", 
                  "706_BOVIS_24HR","716_BOVIS_24HR","721_BOVIS_24HR","724R_BOVIS_25HR","727R_BOVIS_25HR")
  order.rna = c("706cn25","716cn25","721cn25","724Rcn25","727Rcn25","706bovis25","716bovis25","721bovis25","724Rbovis25","727Rbovis25")
  #
  return(
    cor(
      x = as.numeric(rpkm.microarray.intersection[rownames(rpkm.microarray.intersection) == ensembl, order.rna]),
      y = as.numeric(microarray.rnaseq.intersection[microarray.rnaseq.intersection == probeset,order.micro][1,])) # all rows for a probeset have same values, taking 1st
    )
}
# Very long step below!
correlations = apply(X=microarray.rnaseq.intersection, MARGIN=1, FUN=cor.RNA.micro)
save(correlations, file="correlations.RData")
correlations = cbind(microarray.rnaseq.intersection[,c("ID","ensembl_id"], cor=correlations)
save(correlations, file="correlations.RData")

# RNAseq/microarray correlation merged with %probes of probeset aligned
correlations.percAligned = merge(x=correlations, y=modes.percentages, by.x="ID", by.y="probeset", sort=T)
save(correlations.percAligned, file="correlations.percAligned.RData")


# Gets the probesets with the most contradictory data between RNAseq andmicroarray
correlations.percAligned[correlations.percAligned$cor < -0.9,]
############## Worst (=oppposite) RNAseq/microarray correlation with 90% probes aligned
# ID         ensembl_id        cor chromosome      perc total
# Bt.16689.2.A1_at ENSBTAG00000001114 -0.9520225         11 0.9090909    11

plot.RNA.micro = function(probeset, ensembl)
{
  load(file="rpkm.microarray.intersection.RData")
  load(file="microarray.rnaseq.intersection.RData")
  load(file="modes_percentages_total.RData")
  library(bovine.db, quietly=T)
  library(annotate, quietly=T)
  gene.symbol = getSYMBOL(probeset, "bovine.db")
  # Orders the samples for paired correlation
  order.micro = c("706_CN_24HR","716_CN_24HR","721_CN_24HR","724R_CN_25HR","727R_CN_25HR", 
                  "706_BOVIS_24HR","716_BOVIS_24HR","721_BOVIS_24HR","724R_BOVIS_25HR","727R_BOVIS_25HR")
  order.rna = c("706cn25","716cn25","721cn25","724Rcn25","727Rcn25","706bovis25","716bovis25","721bovis25","724Rbovis25","727Rbovis25")
  #
    plot(
      x = as.numeric(rpkm.microarray.intersection[rownames(rpkm.microarray.intersection) == ensembl, order.rna]),
      y = as.numeric(microarray.rnaseq.intersection[microarray.rnaseq.intersection == probeset,order.micro][1,]), # all rows for a probeset have same values, taking 1st
      main =paste(probeset,"=",ensembl,"\n",gene.symbol, "\n",
                  round(modes.percentages[modes.percentages$probeset == probeset,]$perc * 100, digits=2), "% of",
                  modes.percentages[modes.percentages$probeset == probeset,]$total, "probes aligned"),
      xlab="RNAseq: log2(RPKM)",
      ylab="Microarray: log2(intensity)")
}

png(filename="Bt.16689.2.A1_at-ENSBTAG00000001114-PRKD3-inverseCorrelation.png")
plot.RNA.micro(probeset="Bt.16689.2.A1_at", ensembl="ENSBTAG00000001114")
dev.off()

# Gets the probesets with the best correlated data between RNAseq andmicroarray
correlations.percAligned[order(correlations.percAligned$cor, decreasing=T),][1,]
############## Worst (=oppposite) RNAseq/microarray correlation with 90% probes aligned
# ID         ensembl_id        cor chromosome      perc total
# Bt.649.1.S1_at ENSBTAG00000003176 0.998155         29    1    11
png(filename="Bt.649.1.S1_at-ENSBTAG00000003176-JAM3-topCorrelation.png")
plot.RNA.micro(probeset="Bt.649.1.S1_at", ensembl="ENSBTAG00000003176")
dev.off()

# This top correlation shows values over a large range for both microarray and RNAseq
# While the worst correlation further above has large RNAseq variations and small microarray variations
#
# I believe that the correlations should be performed on data normalised between the technologies
# i.e. the variance of RNAseq and microarray should be brought to a comparable level overall

# Plots the distribution of log2(intensity) and log2(RPKM) on the same plot
png(filename="Ranges_RPKM_Intensity.png")
plot(density(rpkm.microarray.intersection[,1]),
     main="Comparison of RNAseq and microarray ranges of expression")
for(i in 2:10){
  lines(density(rpkm.microarray.intersection[,i]))
}
for(i in 2:11){
  lines(density(microarray.rnaseq.intersection[,i]), col="red", lty=2)
}
legend(x="topleft", legend=c("log2(RPKM)", "log2(Intensity)"), col=c("black","red"), lty=c(1,2))
dev.off()
