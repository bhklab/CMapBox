
#' Function to Create a Volcano Plot 
#'
#' This function will create a volcano plot that has the top and bottom 5 genes (according to gene estimates) labelled
#' @param geneEst a numeric vector with the estimates (typically LogFC) for the genes
#' @param genesSymb a character vector with the gene symbols for the gene estimates provided
#' @param genesP a numeric vector with the p values for the gene estimates provided
#' @param idsOfInt a list of Ids to be kept when making the volcano plots, if not provided all Ids will be used
#' @param genesId a character vector with gene Ids for the gene estimates provided. The gene Ids must be of the same type as IdsOfInt. If not provided genes outside of IdsOfInt cannot be removed
#' @keywords gwcCMap
#' @export
#' @examples
#' add example with CMAP here
#' data("geneDataGwc")
#' geneDataClean = cleanData(geneIds = geneDataGwc$ensembl_id, geneEsts = geneDataGwc$logFC, pvals = geneDataGwc$P.Value, forRankAndPlot = TRUE)
#' makeVolcPlot(genesEst = geneDataClean$geneEsts, genesSymb = geneDataClean$symbol, genesP = geneDataClean$pvals)
#' ##only use genes that are present in the pharmacoSet
#' makeVolcPlot(genesEst = geneDataClean$geneEsts, genesSymb = geneDataClean$symbol, genesP = geneDataClean$pvals, genesId = geneDataClean$ensemble, idsOfInt = rownames(drugPertEx))
#' 



makeVolcPlot = function(genesEst, genesSymb, genesP, genesId = NULL, idsOfInt = NULL)
{
  #genesSymb = geneDataDf$symbol
  #genesP = geneDataDf$pvals
  #genesEst = volcPlotEsts
  #IdsOfInt = gsub("_at", "", rownames(drug.perturbation))
  #genesId = gsub("_at", "", geneDataDf$ensemble)
  #idsOfInt = rownames(drugPert)
  #genesId = genesId
  
  if(is.element("calibrate", installed.packages()[,1]) == FALSE)
  {
    install.packages("calibrate")
  }
  library(calibrate)
  
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
  negSubset = volcFrame[mostNegInds[1:5], ]
  mostPosInds = sort(volcFrame[,"geneEst"], index.return = TRUE, decreasing = TRUE)$ix
  posSubset = volcFrame[mostPosInds[1:5], ]

  with(volcFrame, plot(geneEst, -log10(pvalue), pch=20, main="Volcano Plot of Data with Top and Bottom 5 geneEst Genes Labelled", xlim=c(1.2*min(volcFrame$geneEst), 1.2*max(volcFrame$geneEst)), font.axis = 2, cex.lab = 1.3, cex.axis = 1.15))
  with(negSubset, points(geneEst, -log10(pvalue), pch=20, col="blue", cex = 2))
  with(negSubset, textxy(geneEst, -log10(pvalue), labs = geneSymb, cex = 1.5, font = 2, col = "green4"))
  with(posSubset, points(geneEst, -log10(pvalue), pch=20, col="red", cex = 2))
  with(posSubset, textxy(geneEst, -log10(pvalue), labs = geneSymb, cex = 1.5, font = 2, col = "green4"))

}
