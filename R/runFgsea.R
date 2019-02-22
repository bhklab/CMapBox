#' Function to compute a gene set enrichment analysis on a data set
#'
#' This function removes genes from data set one from data set two according to the supplied conditions on the values sent to the function
#' @param pathways list of gene sets to check. If supplying the list ensure the stats variable names are the same as the names within the pathways variable list elements. Enter "all", "bp", "cc", "mf" to have all GO pathways, the biological processes pathways, the cellular component pathways, or the molecular functions pathways used, respectively.
#' @param geneStats the statistics for the genes 
#' @param geneIds the Ids for the genes. Ids must be the same as the ids in the supplied pathways or entrez ids if "all", "bp", "cc", "mf" were used in pathways
#' @param nperm number of permutations to do, default is 50000
#' @return A table with GSEA results. Each row corresponds to a tested pathway. See fgsea for info pertaining to the columns
#' @export
#' @examples
#' data("geneDataGwc")
#' geneDataClean = cleanData(geneIds = geneDataGwc$symbol, geneEsts = geneDataGwc$t, pvals = geneDataGwc$P.Value)
#' geneStatsDat = geneDataClean$geneEsts
#' dataGseaResults = runFgsea(pathways = "cc", geneStats = geneStatsDat, geneIds = geneDataClean$entrez)
#' 
#' @import GSA
#' @import piano
#' @import parallel
#' @import fgsea

runFgsea = function(pathways, geneStats, geneIds, nperm = 50000, numCores = 1)
{
  names(geneStats) = geneIds
  # if(is.element("GSA", installed.packages()[,1]) == FALSE)
  #   install.packages("GSA")
  # library(GSA)
  # if(is.element("piano", installed.packages()[,1]) == FALSE)
  #   install.packages("piano")
  # library(piano)
  # if(is.element("parallel", installed.packages()[,1]) == FALSE)
  #   install.packages("parallel")
  # library(parallel)
  # if(is.element("fgsea", installed.packages()[,1]) == FALSE)
  #   install.packages("fgsea")
  # library(fgsea)
  
  #numCores = 1
  if((detectCores() - 1) > 0)
    coresNum = detectCores() - 1
  
  pathwaysObj = pathways
  if(pathwaysObj[1] == "all")
  {
    data("aGwcPack")
    pathwaysObj = aGwcPack
  }else if(pathwaysObj[1] == "bp"){
    data("abpGwcPack")
    pathwaysObj = abpGwcPack
  }else if(pathwaysObj[1] == "cc"){
    data("accGwcPack")
    pathwaysObj = accGwcPack
  }else if(pathwaysObj[1] == "mf"){
    data("amfGwcPack")
    pathwaysObj = amfGwcPack
  }
  
  gseaRes = fgsea(pathways = pathwaysObj, geneStats, nperm = nperm, nproc = coresNum)
  return(gseaRes)
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

