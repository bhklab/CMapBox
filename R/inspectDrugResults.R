#' Function to inspect drug candidates
#'
#' This function generates plots that will allow one to verify that the drug candidates that should reverse and enhance the disease phenotype are in fact doing so
#' @param genesEst a numeric vector with the estimates (typically LogFC) for the genes
#' @param genesSymb a character vector with the gene symbols for the gene estimates provided
#' @param genesP a numeric vector with the p values for the gene estimates provided
#' @param genesId a character vector with gene Ids for the gene estimates provided
#' @param drugPert the drug perturbation signature used in the analysis
#' @param pharmSet the pharmacoSet used to generate the drug perturbation signature
#' @param drugScoreMeth a string specifying which drug repurposing technique was used to score the drugs during each trial. The options for this parameter are currently "gwc" and "fgsea"
#' @param mDataType the type of molecular data to obtain from the PharmacoSet (see PharmacoGx::molecularProfiles)
#' @param drugResults a drug rank data frame (such as that returned by rankDrugsGwc) ordered by most negative to most positive scored drugs with the rownames containing the drug names
#' @param gwcMethod a character string specifying which method was used when computing correlations in the gwc function. The options are spearman or pearson.
#' @param drugEst a boolean specifying whether the estimates for each gene of the drug perturbation signature were used in the gwc calculation (TRUE) or if the t-stats for each gene in the drug perturbation signature were used in the gwc calculation (FALSE).
#' @param drugVolcPlot a boolean (default is FALSE) specifying whether to generate a volcano plot for the top and bottom 5 genes (by estimate) in the drug that were present in the data supplied
#' @param topDrugName the name of the drug that one would like to check reverses the disease phenotype. If not provided, the best candidate according to the drugResults frame will be selected.
#' @param botDrugName the name of the drug that one would like to check enhances the disease phenotype. If not provided, the best candidate according to the drugResults frame will be selected.
#' @param showMimic a boolean (default is FALSE) specifying whether to show plots for the drug that upregulates and down regulates the genes that are overexpressed and under expressed in the disease state.
#' @export
#' @examples
#' data("geneDataGwc")
#' data("drugPertEx")
#' data("psetSub")
#' data("drugResults")
#' geneDataClean = cleanData(geneIds = geneDataGwc$ensembl_id, geneEsts = geneDataGwc$logFC, pvals = geneDataGwc$P.Value, forRankAndPlot = TRUE)
#' inspectDrugResults(genesEst = geneDataClean$geneEsts, genesSymb = geneDataClean$symbol, genesP = geneDataClean$pvals, genesId = geneDataClean$ensemble, drugPert = drugPertEx, pharmSet = psetSub, mDataType = "rna", drugResults = drugResults, gwcMethod = "pearson", drugEst = TRUE)


inspectDrugResults = function(genesEst, genesSymb, genesP, genesId, drugPert, pharmSet, drugScoreMeth, mDataType, drugResults, gwcMethod, drugEst, drugVolcPlot = FALSE, topDrugName = NULL, botDrugName = NULL, showMimic = FALSE)
{
  # inspectDrugResults(volcPlotEsts, geneDataDf$symbol, geneDataDf$pvals, geneDataDf$ensemble, drugPert, pharmSet, mDataType = mDataType, cmapResults, gwcMethod, drugEst, drugVolcPlot = FALSE, showMimic = showMimic)
  #genesEst = volcPlotEsts[idsUsed]
  #genesSymb = geneDataDf$symbol[idsUsed]
  #genesP = geneDataDf$pvals[idsUsed]
  #genesId = geneDataDf$ensemble[idsUsed]
  #drugResults = cmapResults
  
  makeVolcPlot(genesEst, genesSymb, genesP, genesId, rownames(drugPert))
  
  if(is.null(topDrugName))
    topDrugName = rownames(drugResults)[1]
  driveGenesNegScore = getSigGeneVolcPlot(genesEst, genesSymb, genesP, genesId, drugPert, topDrugName, drugScoreMeth, gwcMethod, droveNegScore = TRUE, drugEst, drugVolcPlot = FALSE, supressPlot = FALSE)
  if(drugScoreMeth == "gwc" | drugScoreMeth == "fgsea")
    print(paste("For the volcano plots that show the genes that drove the negative connectivity score, blue and red labelled genes should cross 0 when going from the data volcano plot to drug", topDrugName,"volcano plot and the blue dots in the volcano plot for the data should be genes that have lower expression in the disease state relative to the normal state. If the prior statements are not true then these are not the drugs you are looking for. You likely need to switch the sign of the supplied geneEsts vector."))
  
  seeDrugImpactGenes(drugPert, pharmSet, mDataType, topDrugName, driveGenesNegScore)
  if(drugScoreMeth == "gwc" | drugScoreMeth == "fgsea")
    print(paste("If the negative gene estimate genes from the data (blue dots in data volcano plots) that drove the negative connectivity score have not had their expression levels increase (Percent Change > 0 and usually blue boxes) on the heatmaps for drug", topDrugName, "then this drug does not appear to increase the expression of genes that had low expression in the disease phenotype relative to the normal phenotype and this is not the drug you are looking for (blue dots in volcano plots should be blue in the heatmap for the drug that reverses the phenotype). "))
    
  if(showMimic == TRUE)
  {
    if(is.null(botDrugName))
      botDrugName = rownames(drugResults)[nrow(drugResults)]
    driveGenesPosScore = getSigGeneVolcPlot(genesEst, genesSymb, genesP, genesId, drugPert, botDrugName, drugScoreMeth, gwcMethod, droveNegScore = FALSE, drugEst, drugVolcPlot = FALSE, supressPlot = FALSE)
    seeDrugImpactGenes(drugPert, pharmSet, mDataType, botDrugName, driveGenesPosScore)
  }
  
  
  
}

