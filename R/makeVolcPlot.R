
#' Function to Create a Volcano Plot
#'
#' This function will create a volcano plot that has the top and bottom 5 genes (according to gene estimates) labelled. Alternatively,
#' users can turn off the plotting option in which case just the info used to make the plot will be returned.
#' @param genesEst a numeric vector with the estimates (typically LogFC) for the genes
#' @param genesSymb a character vector with the gene symbols for the gene estimates provided
#' @param genesP a numeric vector with the p values for the gene estimates provided
#' @param idsOfInt a list of Ids to be kept when making the volcano plots, if not provided all Ids will be used
#' @param genesId a character vector with gene Ids for the gene estimates provided. The gene Ids must be of the same type as IdsOfInt. If not provided genes outside of IdsOfInt cannot be removed
#' @param makePlot a boolean that, if TRUE, will create a volcano plot. Default is TRUE
#' @param topGenes The number of the most positive and negative genes, accordin to geneEst, to label on the volcano plot and/or return. Default is 5
#' @return a list containing the frame used for the volcano plot, a frame with info for the top positive genes labelled in the volcano plot, and
#' a frame with info for the top negative genes labelled in the volcano plot
#' @keywords gwcCMap
#' @export
#' @examples
#' #add example with CMAP here
#' data("geneDataGwc")
#' geneDataClean = mapGenes(geneIds = geneDataGwc$symbol, geneEsts = geneDataGwc$logFC, pvals = geneDataGwc$P.Value, forRankAndPlot = TRUE)
#' volcPlotData = makeVolcPlot(genesEst = geneDataClean$geneEsts, genesSymb = geneDataClean$symbol, genesP = geneDataClean$pvals)
#'



makeVolcPlot = function(genesEst, genesSymb, genesP, genesId = NULL, idsOfInt = NULL, makePlot = TRUE, topGenes = 5)
{
  #genesSymb = geneDataDf$symbol
  #genesP = geneDataDf$pvals
  #genesEst = volcPlotEsts
  #IdsOfInt = gsub("_at", "", rownames(drug.perturbation))
  #genesId = gsub("_at", "", geneDataDf$ensemble)
  #idsOfInt = rownames(drugPert)
  #genesId = genesId

  volcFrame  = as.data.frame(cbind(genesSymb, genesEst, genesP), stringsAsFactors = FALSE)

  if(!is.null(idsOfInt))
  {
    if((is.null(genesId))){
      warning("Could not remove genes that were not of interest because genesId was not provided")
    }else{
      rownames(volcFrame) = genesId
      sameRows = intersect(idsOfInt, genesId)
      volcFrame = volcFrame[row.names(volcFrame) %in% sameRows, ]
    }
  }

  colnames(volcFrame) = c("geneSymb", "geneEst", "pvalue")
  volcFrame = transform(volcFrame, geneEst = as.numeric(geneEst))
  volcFrame = transform(volcFrame, pvalue = as.numeric(pvalue))

  mostNegInds = sort(volcFrame[,"geneEst"], index.return = TRUE)$ix
  negSubset = volcFrame[mostNegInds[1:topGenes], ]
  mostPosInds = sort(volcFrame[,"geneEst"], index.return = TRUE, decreasing = TRUE)$ix
  posSubset = volcFrame[mostPosInds[1:topGenes], ]

  with(volcFrame, plot(geneEst, -log10(pvalue), pch=20, main="Volcano Plot of Data with Top and Bottom 5 geneEst Genes Labelled", xlim=c(1.2*min(volcFrame$geneEst), 1.2*max(volcFrame$geneEst)), font.axis = 2, cex.lab = 1.3, cex.axis = 1.15))
  with(negSubset, points(geneEst, -log10(pvalue), pch=20, col="blue", cex = 2))
  with(negSubset, textxy(geneEst, -log10(pvalue), labs = geneSymb, cex = 1.5, font = 2, col = "green4"))
  with(posSubset, points(geneEst, -log10(pvalue), pch=20, col="red", cex = 2))
  with(posSubset, textxy(geneEst, -log10(pvalue), labs = geneSymb, cex = 1.5, font = 2, col = "green4"))

  retList = list(volcFrame, posSubset, negSubset)
  names(retList) = c("volcFrame", "posSubset", "negSubset")
  return(retList)


}
