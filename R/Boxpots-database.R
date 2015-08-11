## MacroSYS DATABASE 


# Pre-requisites
# A microarray dataset including all probe sets of interest. Here called "farms"
# bioconductor packages, in particular bovine.db and annotate packages

# appends Gene Symbol, EnsEMBL and ENTREZ IDs to the tables of informative probe sets
AnnotateTopTable <- function(df) 
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

##########################
# GeneSymbol <-> ProbeID #
##########################
# Expression values are obtained by probeID
# Database queries are more likely by Gene Symbol

# For the search by Gene Symbol I need a database table converting: GeneSymbol <-> ProbeID

# Select the directory to save the text file in
setwd("C:\\Users\\krue\\Documents\\Kevin-Logs\\20130320_MDM\\Database")

# Retrieve the probe IDs from the experiment of interest
tmp_probes = data.frame(ID=rownames(assayData(farms)$exprs)) # Prepares a dataset with the appropriate number of rows
tmp_probes$ID = rownames(assayData(farms)$exprs) # for AnnotateTopTable() to work, IDs need to be assigned with "tmp_probes$ID = "

# Retrieve the Gene Symbol annotated to each probe ID
tmp_probes_geneSymbols = AnnotateTopTable(tmp_probes)[,c("ID","gene.symbols")]

# Exports the table in a format easy to import in phpmyadmin (field separator '\t', )
write.table(x=tmp_probes_geneSymbols, file="genesymbols_probes.tab", quote=F, sep="\t", eol="\n", dec=".", row.names=F, col.names=F)

# clears memory
rm(tmp_probes, tmp_probes_geneSymbols)

####################
# ProbeID <-> I/NI #
####################
# Normalised expression values are available for all probe sets
# DE calls were only performed on informative (I) probe sets

# For the clarity of the website, I need a database table stating whether a probe set is informative of not

# Select the directory to save the text file in
setwd("C:\\Users\\krue\\Documents\\Kevin-Logs\\20130320_MDM\\MDM5-website\\Database")

# Retrieve the probe IDs from the experiment of interest
tmp_probes = data.frame(ID=rownames(assayData(farms)$exprs)) # Prepares a dataset with the appropriate number of rows

# Labels probe whether they are in the informative data set or not (1 = TRUE, 0 = FALSE)
tmp_probes$I = as.numeric(tmp_probes$ID %in% rownames(assayData(farms_informative)$exprs))

# Exports the table in a format easy to import in phpmyadmin (field separator '\t', )
write.table(x=tmp_probes, file="probes_informative.tab", quote=F, sep="\t", eol="\n", dec=".", row.names=F, col.names=F)

# clears memory
rm(tmp_probes)

####################################
# ProbeID <-> log2FC + adj-p-value #
####################################
# log2FC and adj-p-value are available only for informative probe sets

# To accompany the boxlpots, I need a database table stating whether a probe set is DE ion each contrast or not

# Select the directory to save the text file in
setwd("C:\\Users\\krue\\Documents\\Kevin-Logs\\20130320_MDM\\MDM5-website\\Database")

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

# write it to a file
write.table(x=tmp_probes, file="probes_DEcall.tab", quote=F, sep="\t", eol="\n", dec=".", row.names=F, col.names=F)

# clears memory
rm(tmp_probes)