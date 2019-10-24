
#' Function to map gene ids to entrez IDs, gene symbols and ensemble IDs
#'
#' This function removes genes that have missing geneIDs, geneEsts or pvals, returns the gene symbols, entrez IDs and ensemble IDs for the geneIDs supplied, and removes genes with the same ID by keeping the genes with lower P values (default) or higher abs(estimates)
#' @param geneIds a character vector containing the gene symbols, ensemble IDs, or entrez IDs for the genes to be analyzed
#' @param geneEsts a numeric vector containing estimates for the difference in expression between the two phenotyes. Note that positive values of this estimate (tstat, logFC, etc) must correspond to higher gene expression in the phenotype one would like to reverse
#' @param pvals P values from a t-test assessing the diferential expression of the genes between the two phenotypes
#' @param byP a boolean value specifying whether genes will be kept according to smallest P values (default) or highest abs(estimate)
#' @param forRankAndPlot a boolean value specifying whether the data should be cleaned in a manner that will guarantee it works in the rankDrugsGwc function and volcano plot functions (removes genes that can't be mapped to ensemble IDs)
#' @return a data.frame containing unique gene IDs and their corresponding estimates and p values
#' @export
#' @examples
#' data("geneDataGwc")
#' ##keep genes with the smaller p value of genes with the same id
#' geneDataClean = mapGenes(geneIds = geneDataGwc$symbol, geneEsts = geneDataGwc$logFC, pvals = geneDataGwc$P.Value)
#' ##keep genes with larger abs(logFC) of genes with the same id. Could also supply entrez ids
#' geneDataClean = mapGenes(geneIds = geneDataGwc$symbol, geneEsts = geneDataGwc$logFC, pvals = geneDataGwc$P.Value, byP = FALSE)


mapGenes = function(geneIds, geneEsts, pvals = NULL, byP = TRUE, forRankAndPlot = FALSE)
{

  if(length(geneIds) != length(unique(geneIds)))
    stop("There are duplicated geneIds. Please ensure there is only 1 of each gene in the provided data as we do not know which of the repeat genes should be
          used in the analysis")

  if(is.null(pvals)) pvals = rep(NA, length(geneIds))
  #library(org.Hs.eg.db)

  geneIds = as.character(geneIds)
  geneDataDf = as.data.frame(cbind(geneIds, geneEsts, pvals), stringsAsFactors=FALSE)
  colnames(geneDataDf) = c("geneIds", "geneEsts", "pvals")
  geneDataDf$geneEsts = as.numeric(as.character(geneDataDf$geneEsts))
  geneDataDf$pvals = as.numeric(as.character(geneDataDf$pvals))

  #remove rows with data missing
  notNaRows = rowSums(is.na(geneDataDf)) == 0
  if(sum(notNaRows) < dim(geneDataDf)[1])
  {
    print(paste(dim(geneDataDf)[1] - sum(notNaRows), "genes removed due to NA gene IDs, estimates, and/or p values"))
    geneDataDf =  geneDataDf[notNaRows, ]
  }

  if(length(geneIds) != length(unique(geneIds)))
    stop("There are duplicated geneIds. Please ensure there is only 1 of each gene in the provided data as we do not know which of the repeat genes should be
         used in the analysis")

  #no longer do below, want to make sure user supplies ideal genes for their analysis as oppose to maing decisions for them
  #remove duplicate elements, keep elements with smaller p values as unsure what the supplied estimate is (t-stat, logFC, etc)
  #if(length(unique(geneDataDf$geneIds)) < dim(geneDataDf)[1])
  if(0)
  {

    if(byP == TRUE & sum(is.na(geneDataDf$pvals)) != nrow(geneDataDf))
    {
      #print(paste(dim(geneDataDf)[1] - length(unique(geneDataDf$geneIds)), "genes with the same ID were found. The genes with the smaller p values of those with the same ID will be kept"))
      #keepRows = getUniqueIds(as.character(geneDataDf$geneIds), geneDataDf$pvals, byP)
    }else{
      #print(paste(dim(geneDataDf)[1] - length(unique(geneDataDf$geneIds)), "genes with the same ID were found. The genes with the largest abs(estimates) of those with the same ID will be kept"))
      #keepRows = getUniqueIds(as.character(geneDataDf$geneIds), geneDataDf$geneEsts, byP)
    }
    #geneDataDf = geneDataDf[keepRows, ]
  }

  #
  if(suppressWarnings(sum(is.na(as.numeric(geneIds))) < length(geneIds)/2))
  {
    geneType = "ENTREZID"
    names(geneDataDf)[1] = "entrez"
    mapping = AnnotationDbi::select(org.Hs.eg.db, keytype=geneType, keys=as.character(geneDataDf$entrez), columns=c("ENSEMBL"))
    mapping = mapping[!duplicated(mapping[,1]),]
    geneDataDf$ensemble = mapping[, 2]
    mapping = AnnotationDbi::select(org.Hs.eg.db, keytype=geneType, keys=as.character(geneDataDf$entrez), columns=c("SYMBOL"))
    mapping = mapping[!duplicated(mapping[,1]),]
    geneDataDf$symbol = mapping[, 2]
  } else if(sum(sum(grepl("ENSG00", geneIds))) > length(geneIds)/2){
    geneType  = "ENSEMBL"
    names(geneDataDf)[1] = "ensemble"
    mapping = AnnotationDbi::select(org.Hs.eg.db, keytype=geneType, keys=as.character(geneDataDf$ensemble), columns=c("ENTREZID"))
    mapping = mapping[!duplicated(mapping[,1]),]
    geneDataDf$entrez = mapping[, 2]
    mapping = AnnotationDbi::select(org.Hs.eg.db, keytype=geneType, keys=as.character(geneDataDf$ensemble), columns=c("SYMBOL"))
    mapping = mapping[!duplicated(mapping[,1]),]
    geneDataDf$symbol = mapping[, 2]
  }else{
    geneType = "SYMBOL"
    names(geneDataDf)[1] = "symbol"
    mapping = AnnotationDbi::select(org.Hs.eg.db, keytype=geneType, keys=as.character(geneDataDf$symbol), columns=c("ENTREZID"))
    mapping = mapping[!duplicated(mapping[,1]),]
    geneDataDf$entrez = mapping[, 2]
    mapping = AnnotationDbi::select(org.Hs.eg.db, keytype=geneType, keys=as.character(geneDataDf$symbol), columns=c("ENSEMBL"))
    mapping = mapping[!duplicated(mapping[,1]),]
    geneDataDf$ensemble = mapping[, 2]
  }

  if(forRankAndPlot == TRUE)
  {
    if(sum(is.na(geneDataDf$ensemble)) > 0)
    {
      print(paste(sum(is.na(geneDataDf$ensemble)), "of the unique gene IDs could not be mapped to ensemble IDs and had to be removed"))
      geneDataDf = geneDataDf[-which(is.na(geneDataDf$ensemble)), ]
    }

    if(length(geneDataDf$ensemble) > length(unique(geneDataDf$ensemble)))
    {
      nonUniqueMap = length(geneDataDf$ensemble) - length(unique(geneDataDf$ensemble))
      print(paste("After removing the non unique gene IDs and unique gene IDs that did not map to ensemble IDs", nonUniqueMap, "of the unqiue gene IDs mapped to the same ensemble ID and were removed using the p values/estimates supplied"))
      if(byP == TRUE)
      {
        keepRows = getUniqueIds(as.character(geneDataDf$ensemble), geneDataDf$pvals, byP)
      }else{
        keepRows = getUniqueIds(as.character(geneDataDf$ensemble), geneDataDf$geneEsts, byP)
      }
      geneDataDf = geneDataDf[keepRows, ]
    }
  }

  return(geneDataDf)

}
