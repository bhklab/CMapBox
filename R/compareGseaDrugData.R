#' Function to determine pathways that are the most significantly enriched in both a data set and in a drug present in a pharmaco set
#'
#' This function computes gsea results for both a data set and a drug. The function returns a table of the results for the pathways that is ordered according to the most significant pathways in both the data and the drug.
#' Significance is determined by ranking the pathways p adjusted for the drug and the data (lower p adjusted --> lower rank) and then averaging the rank of the pathway in both the drug and the data. Lower average ranks appear higher in the data frame returned from this function.
#' @param pathways list of gene sets to check. Gene IDs in the gene sets must be entrez IDs or ensemble IDs
#' @param geneStats the statistics for the genes with names the same as the ids in the supplied pathways or entrez ids if "all", "bp", "cc", "mf" were used in pathways
#' @param drugName the name of the drug in the pharmacoSet to perform the analysis on
#' @param drugPert a drug perturbation signature (object of class PharmacoSig) with rownames that correspond to either entrez IDs or ensemble IDs for the genes.
#' @param useTstat a boolean specifying whether to use the drugs t-stat in fgsea (default is TRUE). If FALSE the drug's estimate for each gene will be used
#' @param nperm number of permutations to do, default is 50000
#' @return A table with GSEA results for the drug and data ordered from most significant pathways in the data and the drug to least. Each row corresponds to a tested pathway. See fgsea for info pertaining to the columns
#' @export
#' @examples
#' data("drugPertEx")
#' data("geneDataGwc")
#' geneDataClean = cleanData(geneIds = geneDataGwc$symbol, geneEsts = geneDataGwc$t, pvals = geneDataGwc$P.Value)
#' geneStatsDat = geneDataClean$geneEsts
#' pathOverlapTab = compareGseaDrugData(pathways = "cc", geneStats = geneStatsDat, geneIds = geneDataClean$entrez, drugName = "BRD-K78431006", drugPert = drugPertEx)

compareGseaDrugData = function(pathways, geneStats, geneIds, drugName, drugPert, useTstat = TRUE, nperm = 50000)
{
  gseaResData = runFgsea(pathways, geneStats, geneIds, nperm = nperm)
  gseaResDrug = runFgseaPharmDrug(pathways, drugName, drugPert, useTstat = useTstat, nperm = nperm)
  
  overlapPaths = intersect(gseaResData$pathway, gseaResDrug$pathway)
  dataIndsKeep = which(gseaResData$pathway %in% overlapPaths)
  gseaResData = gseaResData[dataIndsKeep, ]
  drugIndsKeep = which(gseaResDrug$pathway %in% overlapPaths)
  gseaResDrug = gseaResDrug[drugIndsKeep, ]
  
  dataPathRanks = rank(gseaResData$padj, ties.method = "min")
  drugPathRanks = rank(gseaResDrug$padj, ties.method = "min")
  avgPathRank = (dataPathRanks + drugPathRanks)/2
  rowOrd = order(avgPathRank, decreasing = FALSE)
  gseaResData = gseaResData[rowOrd, ]
  gseaResDrug = gseaResDrug[rowOrd, ]
  colnames(gseaResData) = paste0("data_", colnames(gseaResData))
  colnames(gseaResDrug) = paste0("drug_", colnames(gseaResDrug))
  compFrame = cbind(gseaResData$data_pathway, sort(avgPathRank, decreasing = FALSE))
  colnames(compFrame) = c("pathway name", "average rank of p adjusted")
  gseaResData = as.data.frame(gseaResData)
  gseaResDrug = as.data.frame(gseaResDrug)
  for(i in 2:ncol(gseaResData))
  {
    compFrame = cbind(compFrame, gseaResData[, as.numeric(i)])
    colnames(compFrame)[ncol(compFrame)] = colnames(gseaResData)[as.numeric(i)]
    compFrame = cbind(compFrame, gseaResDrug[, as.numeric(i)])
    colnames(compFrame)[ncol(compFrame)] = colnames(gseaResDrug)[as.numeric(i)]
  }
  
  return(compFrame)
}

#setwd("C:\\Users\\micha\\Documents\\PMH Research\\Code from Neel Project")
#geneData = read.table("breast_normal_vs_dtp-all_cell_lines.txt",sep="\t", header=TRUE)
#geneIds = as.character(geneData$ensembl_id)
#geneEsts = as.numeric(as.character(geneData$t))
#pvals = as.numeric(as.character(geneData$P.Value))
#geneDataDf = cleanData(geneIds, geneEsts, pvals, forDrugRank = TRUE)
#geneStats = geneDataDf$geneEsts
#names(geneStats) = geneDataDf$entrez
#gseaResData = runFgsea("cc", geneStats)
#setwd("C:\\Users\\micha\\Documents\\PMH Research\\GWC CMAP Package")
#load("cmap_sig_rna.RData")
#gseaResDrug = runFgseaPharmDrug("cc", "vinburnine", drug.perturbation, useTstat = TRUE, nperm = 50000)
