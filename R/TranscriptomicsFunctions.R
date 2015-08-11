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
### Test data frames ###
#test = data.frame(ID=c("Bt.17794.1.S1_at","Bt.523.1.S1_at"), stringsAsFactors=FALSE)
#fail1 = data.frame(IDe=c("Bt.17794.1.S1_at","Bt.523.1.S1_at"), stringsAsFactors=FALSE)
### Test usage ###
#AnnotateProbesTable(df=test)
#AnnotateProbesTable(df=test, ENSEMBL=F)
#AnnotateProbesTable(df=test, ENTREZ=F)
#AnnotateProbesTable(df=test, gene.symbol=F,  ENSEMBL=T, ENTREZ=F)
#AnnotateProbesTable(df=fail1, probeCol="IDe")


EntrezEnsemblStat = function(df, annotPkg="bovine.db", probeCol="ID")
  {
  # Counts the number of probes havingENTREZ or ENSEMBL annotations in a dataset.
  #
  # Args:
  #   df: Data frame containing the probe ID.
  #   annotPkg: Annotation package to use. Default is "bovine.db"
  #   probeCol: identifier of the column containing the probe IDs. Default is "ID".
  #
  # Returns:
  #   A table with counts and percentage of annotated genes.
  #
  total = nrow(df)
  # ENSEMBL annotation statistics
  tmp = AnnotateProbesTable(df, annotPkg=annotPkg, probeCol=probeCol,
                   gene.symbol=FALSE, ENSEMBL=TRUE, ENTREZ=FALSE)
  Ensembl.count = sum(as.numeric(
    lapply(
      aggregate(x=tmp, by=list(tmp[,probeCol]),
                FUN=function(x) !is.na(x))$ensembl_id, FUN=any)
    ))
  # ENTREZ annotation statistics
  tmp = AnnotateProbesTable(df, annotPkg=annotPkg, probeCol=probeCol,
                         gene.symbol=FALSE, ENSEMBL=FALSE, ENTREZ=TRUE)
  Entrez.count = sum(as.numeric(
    lapply(
      aggregate(x=tmp, by=list(tmp[,probeCol]),
                FUN=function(x) !is.na(x))$gene_id, FUN=any)
  ))
  return(data.frame(cbind(Count=c(Entrez.count, Ensembl.count),
                   Percentage=c(Entrez.count/length(df$ID), Ensembl.count/length(df$ID)),
                   Total=rep(length(df$ID),2)),
             row.names=c("Entrez","Ensembl")))
}
### Test data frames ###
#test = data.frame(ID=c("Bt.17794.1.S1_at","Bt.523.1.S1_at"), stringsAsFactors=FALSE)
#fail1 = data.frame(IDe=c("Bt.17794.1.S1_at","Bt.523.1.S1_at"), stringsAsFactors=FALSE)
### Test usage ###
#EntrezEnsemblStat(df=test)
#EntrezEnsemblStat(df=fail1, probeCol="IDe")

UpsAndDowns = function(df, pCol="adj.P.Val", logFCcol="logFC", pThresh=0.05)
{
  # Counts the number of probes having ENTREZ or ENSEMBL annotations in a dataset.
  #
  # Args:
  #   df: Data frame containing the differential expression data.
  #   pCol: identifier of the column containing the (adjusted) p-value. Default is "adj.P.Val".
  #   logFCcol:  identifier of the column containing the log fold-change. Default is "logFC".
  #   pThresh:  Threshold on the (adjusted) p-value. Default is 0.05.
  #
  # Returns:
  #   A table with counts of probesets up and down-regulated.
  #
  T = nrow(df)
  # Error handling
  if (!logFCcol %in% colnames(df)){
    stop("Column does not exist in data.frame: ",
         logFCcol, ".")
  }
  if (!pCol %in% colnames(df)){
    stop("Column does not exist in data.frame: ",
         pCol, ".")
  }
  if (pThresh <= 0 | pThresh > 1){
    stop("Invalid value for threshold on (adjusted) p-value: ",
         pThresh, ".")
  }
  # Count number of probe sets significantly up-regulated
  u = nrow(df[df$adj.P.Val < pThresh & df$logFC > 0,])
  # Count number of probe sets significantly down-regulated
  d = nrow(df[df$adj.P.Val < pThresh & df$logFC < 0,])
  # Total of differentially expressed genes
  t = u + d
  # Count number not significant
  n = nrow(df[df$adj.P.Val >= pThresh,])
  # Puts the results together
  res = data.frame(cbind(Ups=u, Downs=d, total.DE=t, not.signf=n, Total=T), row.names="Counts")
  return (res)
}
### Test data frames ###
#test = data.frame(adj.P.Val=c(0.01,0.06), logFC=c(-0.1, 0.5), stringsAsFactors=FALSE)
#fail1 = data.frame(IDe=c("Bt.17794.1.S1_at","Bt.523.1.S1_at"), stringsAsFactors=FALSE)
### Test usage ###
#UpsAndDowns(df=test)
#UpsAndDowns(df=test, pThresh=0.1)


WriteContrastsFromBayes = function (eB, cntrs,
                                      annotPkg = "bovine.db", gene.symbol=TRUE,
                                      ENTREZ=FALSE, ENSEMBL=TRUE,
                                      number=NA, adj.method="BH", sort.by="P",
                                      quote=FALSE, sep="\t", eol="\n", na="NA",
                                      row.names=FALSE, col.names=TRUE,
                                      out.dir=NA)
{
  # Extract the top-ranked genes from a linear model fit for each of the 
  # specified contrasts, and writes them to respective files.
  #
  # Args:
  #   eB: eBayes variable containing the differential expression data.
  #   cntrs: List of contrasts to search in the dataset.
  #   gene.symbol: If TRUE, annotates output files with official gene symbol.
  #   ENSEMBL: If TRUE, annotates output files with ENSEMBL gene ID.
  #   ENTREZ: If TRUE, annotates output files with ENTREZ gene ID.
  #   number: Maximal number of rows in the files produced. (NOTE: Excluding header)
  #             Default is the number of rows in the eBayes dataset.
  #   adj.method: Method used to computed an adjusted p-value. Default is "BH".
  #   sort.by:  Character string specifying statistic to rank genes by.
  #   quote:  Whether to quote values in output files. Default is "FALSE".
  #   sep:  The field separator string. Default is "\t".
  #   eol: The character(s) to print at the end of each line (row). Default is "\n".
  #   na: The string to use for missing values in the data. Default is "NA".
  #   row.names:  Whether to write row names in output files. Default is "FALSE".
  #   col.names:  Whether to write column names in output files. Default is "TRUE".
  #   out.dir:  Directory where to store the output files.
  #   annotPkg: Annotation package to use. Default is "bovine.db"
  
  #   
  # Returns:
  #   NA
  #
  n = nrow(eB)
  #Error handling
  if (length(cntrs) == 0){
    stop("No element in list of contrast names: cntrs")
  }
  if (!is.na(number) & number < 0){
    stop("Invalid value for number of rows in output file: ",
         number, ".")
  }
  if (!adj.method %in% c("none", "BH", "BY", "holm")){
    stop("Invalid value for p-value adjustment method: ",
         adj.method, ".")
  }
  if (!sort.by %in% c("logFC", "AveExpr", "t", "P", "p", "B", "none", 
                      "M","A","Amean","T","p")){
    stop("Invalid value for statistic to sort by: ",
         sort.by, ".")
  }
  # If any, returns the first contrast name which cannot be found in the dataset
  for (cntr in cntrs){
    if (!cntr %in% colnames(eB$contrasts)){
      stop("Contrast name cannot be found in dataset: ",
           cntr, ".")
    }
  }
  # If not user-defined, sets the "number" value
  if (is.na(number)){
    number = n
  }
  # If not defined, sets output directory to current one
  if (is.na(out.dir)){
    out.dir = getwd()
  }
  # For each contrast, produces the files
  for(cntr in cntrs){
    DE = topTable(eB, coef= cntr, number=number, adjust.method=adj.method, sort.by=sort.by)
    if (any(gene.symbol, ENSEMBL, ENTREZ)){
      DE = AnnotateProbesTable(DE, annotPkg=annotPkg, gene.symbol=gene.symbol, ENSEMBL=ENSEMBL, ENTREZ=ENTREZ)
    }
    write.table(DE, file = paste(out.dir,"/",cntr,".txt", sep=""), quote = quote, sep = sep, eol = eol, na = "NA", row.names = row.names, col.names = col.names)
    # Prints a message for each file written
    cat("Wrote file:", paste(out.dir,"/",cntr,".txt", sep=""), fill=TRUE)
  }
}
### Test data frames ###
# To complicated to create fake one.
### Test usage ###
#ExtractContrastsFromBayes(eB=eb_filtered, cntrs=list("BOVIS_2HR-MAP_2HR","BOVIS_6HR-MAP_6HR"),
#                          annotPkg="bovine.db", gene.symbol=TRUE, ENTREZ=FALSE, ENSEMBL=FALSE,
#                          number=NA, adj.method="BH", sort.by="P", quote=FALSE, sep="\t", eol="\n",
#                          na="NA", row.names=FALSE, col.names=TRUE, out.dir="C:/Users/krue/Desktop/test")

DEthesholdsCounts = function(eB, cntrs,
                             thresh=list(FCs=c(0), pVals=c(0.05)))
{
  # Counts the number of probes satisfying given threshold for differential expression
  # in given contrasts in a linear model dataset.
  #
  # Args:
  #   eB: eBayes variable containing the differential expression data.
  #   cntrs: List of contrasts to search in the dataset.
  #   thresh: List of fold-changes and p-values to use as threshold for differential
  #           expression.  
  #
  # Returns:
  #   A table with counts of probesets up, down-regulated, total differentially 
  #   expressed, not-differentially expressed, and total in dataset.
  #
  n = length(thresh$FCs)
  N = nrow(eB)
  #Error handling
  if (length(cntrs) == 0){
    stop("No element in list of contrast names: cntrs")
  }
  if (length(thresh$FCs) != length(thresh$pVals)){
    stop("FCs and pVals element of thresh argument are of different length: ",
         length(thresh$FCs), " and ", length(thresh$pVals), ".")
  }
  if (sum(thresh$pVals <= 0) > 0){
    stop("At least one invalid p-value: <= 0")
  }
  if (sum(thresh$pVals > 1) > 0){
    stop("At least one invalid p-value: > 1")
  }
  # If any, returns the first contrast name which cannot be found in the dataset
  for (cntr in cntrs){
    if (!cntr %in% colnames(eB$contrasts)){
      stop("Contrast name cannot be found in dataset: ",
           cntr, ".")
    }
  }
  # Initialises the column names of the output table with the total number of probesets
  cNames = c("Total")
  for (i in 1:n){
    substr1 = paste(as.character(thresh$FCs[i]), "_",as.character(thresh$pVals[i]), sep="")
    cNames = append(x=cNames, values=c(
      paste(substr1, "_up", sep=""),
      paste(substr1, "_down", sep=""),
      paste(substr1, "_DE", sep=""),
      paste(substr1, "_not.sign", sep="")))
  }
  # For each contrast, produces the corrresponding file
  DE = data.frame()
  for(cntr in cntrs){
    tT = topTable(fit=eB, coef= cntr, number=N)
    # Initisalises the list of counts with the total number of probesets
    counts = c(N)
    for (i in 1:n){
      # Get the Up and Down information for each threshold for a given contrast
      UnD = UpsAndDowns(df=tT, pCol="adj.P.Val", logFCcol="logFC", pThresh=0.05)
      # Append data for each threshold in the array for a given contrast
      counts = append(x=counts, values=c(
        UnD$Ups, UnD$Downs, UnD$total.DE, UnD$not.signf))
    }
    DE = rbind(DE, counts)
  }
  colnames(DE) = cNames
  rownames(DE) = cntrs
  return(DE)
}
### Test data frames ###
# To complicated to create fake one.
### Test usage ###
#DEthesholdsCounts(eB=eb_filtered, cntrs=list("BOVIS_2HR-MAP_2HR","BOVIS_6HR-MAP_6HR"),
#                  thresh=list(FCs=c(0,0.1), pVals=c(0.05,0.01)))

are.DEprobes.contrasts = function(eB, cntrs=colnames(eB), probes=rownames(eB), pThresh=0.05,
                                  adj.method="BH")
{
  # Returns whether probesets are DE in given contrasts.
  #
  # Args:
  #   eB: eBayes variable containing the differential expression data.
  #   cntrs: List of contrasts to search in the dataset. Defaults to all.
  #   probes: List of probesets to search in the dataset. Defaults to all.
  #   pThresh: Threshold on adjusted p-value for 
  #           expression. Default is 0.05.
  #   adj.method: Method used to adjust the p-values for multiple testing.
  #               Any of "none", "BH", "BY" and "holm". Default is "BH".
  #
  # Returns:
  #   A table of TRUE/FALSE whether given probes are DE in givent contrasts
  #
  N = nrow(eB)
  #Error handling
  if (length(cntrs) == 0){
    stop("No element in list of contrast names: ", cntrs)
  }
  if (length(probes) == 0){
    stop("No element in list of probesets: ", probes)
  }
  if (pThresh <= 0 | pThresh > 1){
    stop("Invalid p-value threshold: ", pThresh)
  }
  if (!adj.method %in% list("none","BH","BY","holm")){
    stop("No element in list of contrast names: cntrs")
  }
  # Prepares the output table
  result = data.frame(ID=probes, stringsAsFactors=FALSE)
  # For each contrast, returns the status of the genes of interest
  for (cntr in cntrs){
    # Obtains the information for all probes for a given contrast
    topT = topTable(eb_filtered, coef= cntr, number=N, adjust.method=adj.method)
    # Filters for required probesets
    topT = topT[topT$ID %in% probes,]
    # Boolean whether the genes are DE or not
    topT[,cntr] = topT$adj.P.Val < 0.05
    # Merges the boolean with the output table
    result = merge(x=result, y=topT[,c("ID",cntr)], by.x="ID", by.y="ID")
  }
  rownames(result) = result$ID
  return(result)
}
### Test data frames ###
# To complicated to create fake one.
### Test usage ##
#are.DEprobes.contrasts(eB=eb_filtered, cntrs=c("BOVIS_2HR-MAP_2HR","BOVIS_6HR-MAP_6HR"),
#                       probes=c("AFFX-BioB-3_at","Bt.188.1.S1_at"), pThresh=0.05, adj.method="BH")

boxplot_expression = function(eSet, probeset, groupCol = "Group", groups="",
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
          cex.axis=0.75, las=3,
          main=paste(probeset, "=", gene.symbol, "| Expression plot"),
          col=col)
  # If requested return the values plotted
  if(values){
    return(df)
  }
}
### Test data frames ###
# To complicated to create fake one.
### Test usage ###
# boxplot_expression(eSet=farms, probeset="Bt.188.1.S1_at",
#                groups=c("CN_0HR","CN_2HR","CN_6HR","CN_24HR","MAP_2HR","MAP_6HR","MAP_24HR",
#                         "BCG_2HR","BCG_6HR","BCG_24HR","BOVIS_2HR","BOVIS_6HR","BOVIS_24HR"),
#                annotPkg="bovine.db", ylim=c(5,15), cex.axis=0.78, las=3,
#                col=c(rep("grey",4),rep("lightblue",3),rep("green",3),rep("red",3)), values=TRUE)
# 
# boxplot_expression(eSet=farms.informative.animal.corrected, probeset="Bt.188.1.S1_at", groupCol="Group",
#                 groups= c("CN_2HR","CN_6HR","CN_24HR","MAP_2HR","MAP_6HR","MAP_24HR","BCG_2HR","BCG_6HR","BCG_24HR",
#                           "BOVIS_2HR","BOVIS_6HR","BOVIS_24HR"),
#                 annotPkg="bovine.db", ylim=c(5,15), cex.axis=0.78, las=3,
#                 col=c(rep("grey",3),rep("lightblue",3),rep("green",3),rep("red",3)), values=FALSE)

ExtractContrastsFeatureFromBayes = function (eB, cntrs=colnames(eB), column="logFC",
                                      adj.method="BH", probesets=rownames(eB))
{
  # Extract a metric from each of a list of specified contrasts.
  #
  # Args:
  #   eB: eBayes variable containing the differential expression data.
  #   cntrs: List of contrasts to search in the dataset.
  #   column: Column to extract in each contrast result table
  #   adj.method: Method used to computed an adjusted p-value. Default is "BH"
  #   probesets:  Subset of probesets to return data about.  
  #   
  # Returns:
  #     An aggregated data frame.
  #
  n = nrow(eB)
  #Error handling
  if (length(cntrs) == 0){
    stop("No element in list of contrast names: cntrs")
  }
  if (!adj.method %in% c("none", "BH", "BY", "holm")){
    stop("Invalid value for p-value adjustment method: ",
         adj.method, ".")
  }
  if (!column %in% c("ID","logFC","AveExpr","t","P.Value","adj.P.Val","B")){
    stop("Invalid value for column to extract: ",
         col, ".")
  }
  # If any, returns the first contrast name which cannot be found in the dataset
  for (cntr in cntrs){
    if (!cntr %in% colnames(eB$contrasts)){
      stop("Contrast name cannot be found in dataset: ",
           cntr, ".")
    }
  }
  #Initialises the output data frame
  result = data.frame(row.names=probesets)
  result$ID = probesets
  # For each contrast, obtains the requested column
  for(cntr in cntrs){
    # Extracts the table of stastistics from the linear model
    DE = topTable(fit=eB, coef= cntr, adjust.method=adj.method, number=n)
    # Renames the column to extract with the name of the contrast to make it unique in the output table
    names(DE)[names(DE) == column] = cntr
    # Merges the column to the output data frame
    result = merge(x=result, y=DE[,c("ID",cntr)], by="ID", sort="T", all.x=TRUE)
  }
  result$ID = NULL
  rownames(result) = probesets
  return(result)
}

### Test data frames ###
# To complicated to create fake one.
### Test usage ###
#plot(ExtractContrastsFeatureFromBayes(eB=eb_filtered, cntrs=c("BOVIS_2HR-CN_2HR", "MAP_2HR-CN_2HR","BCG_2HR-CN_2HR" )))

SharedDEprobesets = function(eB, direction, cntrs=colnames(eB), logFC=0, p.value=0.05, adj.method="BH")
{
  # Filters for probesets with differential expression in given 
  #
  # Args:
  #   eB: eBayes variable containing the differential expression data.
  #   direction: One of "up" "down"-regulation. Should the genes be shared up- or down-regulated ones?
  #   cntrs: List of contrasts to search in the dataset.
  #   logFC:  Shared threshold on log fold change.
  #   p.value:  Shared threshold on p-value.
  #   adj.method: Method used to computed an adjusted p-value. Default is "BH".
  
  #   
  # Returns:
  #     An aggregated data frame, including probe ID, fold changes and p-values.
  #  
  n = nrow(eB)
  # Error handling
  if (!direction %in% c("up","down")){
    stop("Invalid value for direction of regulation: ", direction)
  }
  if (length(cntrs) == 0 | any(cntrs == "")){
    stop("No element in list of contrast names: cntrs")
  }
  if (p.value <= 0 | p.value > 1){
    stop("Invalid p-value threshold: ", pThresh)
  }
  if (!adj.method %in% c("none", "BH", "BY", "holm")){
    stop("Invalid value for p-value adjustment method: ",
         adj.method, ".")
  }
  # Initialises the output data frame
  result = data.frame(row.names=rownames(eB))
  result$ID = rownames(result)
  for(cntr in cntrs){
    topT = topTable(eB, coef= cntr, number=n, adjust.method=adj.method)[,c("ID","logFC", "adj.P.Val")]
    if(direction == "up"){
      topT = topT[topT$logFC > logFC & topT$adj.P.Val < p.value,]
    }
    else{ # we already ensured that direction is either "up" or "down"
      topT = topT[topT$loFC < logFC & topT$adj.P.Val < p.value,]
    }
    colnames(topT)[colnames(topT)=="logFC"] = paste("logFC.", cntr, sep="")
    colnames(topT)[colnames(topT)=="adj.P.Val"] <- paste("adj.P.Val.", cntr, sep="")
    result = merge(x=result, y=topT, by="ID", sort=TRUE)
  }
  rownames(result) = result$ID
  return(result)
}