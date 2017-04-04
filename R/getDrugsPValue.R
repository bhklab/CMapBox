
#' Function to compute p values for the drugs
#'
#' This function will compute the p values for only the drugs inputted via drugNameVec. To save time, users are recommended to run rankDrugsGwc with a low nperm and then use this function with a high nperm (~100,000) to get the p values for the top drug candidates.
#' @param drugNameVec a character vector containing the names of the drugs in the drug perturbation signature to compute a more accurate p value for
#' @param geneIds a character vector containing the gene symbols, ensemble IDs, or entrez IDs for the genes to be analyzed
#' @param geneEsts a numeric vector containing estimates for the difference in expression between the two phenotyes. Note that positive values of this estimate (tstat, logFC, etc) must correspond to higher gene expression in the phenotype one would like to reverse
#' @param pvals P values from a t-test assessing the diferential expression of the genes between the two phenotypes
#' @param drugPert a drug perturbation signature (object of class PharmacoSig) with rownames that correspond to the ensemble IDs of the genes
#' @param drugScoreMeth a string specifying which drug repurposing technique to use to score the drugs during each trial. The options for this parameter are currently "gwc" and "fgsea". Default is "gwc"
#' @param pCut a boolean specifying whether to use pvalues to remove insignificant genes from the CMAP analysis (TRUE) or to remove genes from the CMAP analysis according to their supplied gene estimates (FALSE)
#' @param cutOff if pCut is TRUE then this value represents the p value threshold used to filter out genes. If pCut is FALSE then this value represents the fraction of genes present in both the data and drug perturbation signature with the top absolute value of gene estimates that will be left in the analysis. cutOff should be between 0 and 1.
#' @param genesToStart a value between 0 and 1 representing the fraction of significant genes present in the data and drug perturbation signature to use in the first iteration of the CMAP analysis. Recommended to be at least 0.15 to avoid large changes during the early iterations due to having too few genes present.
#' @param numbIters The number of iterations that will occur in the analysis. numbIters will set the rate at which genes are added to the analysis.
#' @param numbPerms The number of permutations to be used to compute the p value in the gwc function. 
#' @param gwcMethod a character string specifying which method to use when computing correlations in the gwc function. The options are spearman (default) or pearson.
#' @param drugEst a boolean specifying whether to use the estimates for each gene of the drug perturbation signature in the gwc calculation (TRUE) or to use the t-stats for each gene in the drug perturbation signatur ein the gwc calculation (FALSE). Default is TRUE
#' @param extraData a data.frame where each column represents values one would like to inspect to determine if the genes corresponding to these values should be removed if they meet the conditions specified in extraCut and extraDirec. Useful if one would like to remove genes based on logFC or other values. extra data must have the same number of rows as the length of vectors geneIds, geneEsts, and pvals.
#' @param extraCut a numeric vector, with the nth value corresponding to the nth column of the m by n data.frame supplied in extraData, where each value indicates the value that needs to be reached for the columns of the data fram in order for the ids with that value to be kept or removed. Whether removal occurs when a value is greater or less than the value in this vector is specified in the extraDirec variable.
#' @param extraDirec a boolean vector, where TRUE means that ids whose value is greater than that specified in extraVals will be removed and FALSE means those with a value less than conditionVals will be removed. 
#' @return a data frame containing the p value from the connectivity score computed using the number of genes specified in the rows for each drug. 
#' @keywords gwcCMap
#' @export
#' @examples
#' data("geneDataGwc")
#' data("drugPertEx")
#' geneDataClean = cleanData(geneIds = geneDataGwc$symbol, geneEsts = geneDataGwc$t, pvals = geneDataGwc$P.Value)
#' pFrame = getDrugsPValue(drugNameVec = c("BRD-K78431006", "BRD-A68739437", "BRD-K73368362"), geneIds = geneDataClean$ensemble, geneEsts = geneDataClean$geneEsts, pvals = geneDataClean$pvals, drugPert = drugPertEx, numbIters = 5, numbPerms = 50000)

getDrugsPValue = function(drugNameVec, geneIds, geneEsts, pvals, drugPert = NULL, drugScoreMeth = "gwc", pCut = TRUE, cutOff = 0.05, genesToStart = 0.20, numbIters = 10, gwcMethod = "spearman", numbPerms = 100000, drugEst = TRUE, extraData = NULL, extraCut = NULL, extraDirec = NULL)
{
  #convert geneEsts and pvals in case they are factors or character vectors
  pharmSet = NULL
  geneEsts = as.numeric(as.character(geneEsts))
  pvals = as.numeric(as.character(pvals))
  
  geneDataDf = cleanData(geneIds, geneEsts, pvals, forRankAndPlot = TRUE)
  
  #pharmacoGx doesnt appear to install if needed, below installs it if it isn't been before
  if(is.element("PharmacoGx", installed.packages()[,1]) == FALSE)
  {
    source("https://bioconductor.org/biocLite.R")
    biocLite("PharmacoGx")
  }
  library(PharmacoGx)
  if(is.null(drugPert))
  {
    print("No drug data was supplied. The drug data present in CMAP will be downloaded and used in the analysis")
    #should just download perturbation signature directly if its available on the website or supply it with the package
    PharmacoGx::downloadPSet("CMAP")
    print("Computing perturbation signature from CMAP data.....this may take half a day, considering finding and supplying a precomputed pset")
    drugPert = drugPerturbationSig(CMAP, mDataType="rna")
    print("Perturbation signature computation finished")
  }
  
  #deal with cmap rownames
  if(grepl("_at", rownames(drugPert)[1]))
    rownames(drugPert) = gsub("_at_", "", paste0(rownames(drugPert), "_"))
  #geneDataDf$ensemble = paste0(geneDataDf$ensemble, "_at")
  
  #deal with L1000 dataset that is ordered different from cmap data and has entrez ids instead of ensemble ids
  if(dim(drugPert)[3] > dim(drugPert)[2])
    drugPert = aperm(drugPert, c(1,3,2))
  if(grepl("geneid.", rownames(drugPert)[1]) == TRUE)
  {
    geneEntrezs = gsub("geneid.", "", rownames(drugPert))
    mapping = AnnotationDbi::select(org.Hs.eg.db, keytype="ENTREZID", keys=as.character(geneEntrezs), columns=c("ENSEMBL"))
    mapping = mapping[!duplicated(mapping[,1]),]
    rownames(drugPert) = mapping[, 2]
  }
  
  #remove genes that are not present in both the data and drugPert
  presentVec = c()
  for(i in 1:nrow(drugPert))
    presentVec = c(presentVec, which(geneDataDf$ensemble == rownames(drugPert)[i]))
  print(paste(length(presentVec),"of the", nrow(geneDataDf),"genes from the data are present in the drug perturbation signature and will be used in the analysis"))
  geneDataDf = geneDataDf[presentVec, ]
  
  genesCmapAll = as.data.frame(NULL)
  genesCmapAll[1:nrow(geneDataDf), 1] = as.numeric(as.character(geneDataDf$geneEsts))
  genesCmapAll[1:nrow(geneDataDf), 2] = as.numeric(as.character(geneDataDf$pvals))
  rownames(genesCmapAll) = geneDataDf$ensemble
  colnames(genesCmapAll) = c("estimate", "pvalue")
  
  if(cutOff > 1)
  {
    while(cutOff > 1)
      cutOff = cutOff/100
    
    warning("The cutOff value supplied was greater than 1 and has been reduced to a fraction between 0 and 1")
    if(pCut == TRUE)
      print(paste("If it was intended that genes with a p value less than", cutOff, "should be kept, then the analysis will be conducted as intended"))
    if(pCut == FALSE)
      print(paste("If it was intended that genes with an abs(estimate) in the top", cutOff*100, " percent should be kept, then the analysis will be conducted as intended"))
  }
  
  if(pCut == TRUE)
  {
    sigGenes = which(genesCmapAll$pvalue < cutOff)
    print(paste(length(sigGenes), "significant genes/genes with p <", cutOff,"in the analysis"))
  }else{
    estOrd = order(abs(genesCmapAll$estimate), decreasing = TRUE)
    numKeep = floor(length(estOrd)*cutOff)
    sigGenes = fcOrd[1:numKeep]
    print(paste(cutOff,"percent of genes with the top esimate will be used in the analysis for a total of", length(sigGenes), "genes"))
  }
  
  genesCmap = genesCmapAll[sigGenes,]
  
  if(!is.null(extraData))
  {
    if(is.null(dim(extraData)))
    {
      extraData = cbind(extraData)
    }
    extraData = extraData[as.numeric(rownames(geneDataDf))]
    rowsToKeep = removeGenes(idsOne = rownames(genesCmapAll), idsTwo = rownames(genesCmapAll), valuesTwo = extraData, conditionVals = extraCut, conditionsDirec = extraDirec)
    if(length(rowsToKeep) > 0){
      genesCmap = genesCmap[rowsToKeep, ]
      print(paste("extraData was supplied and used to remove an additional ",length(sigGenes) - length(rowsToKeep) ," genes, the number of genes left in the analysis now is ", length(rowsToKeep)))
    }else
    {
      warning("no genes were left after using the criteria specified by extraData, extraCut, and extraDirec so the analysis will proceed without using the extra criteria")
    }
  }
  
  #when adding genes during analysis, add them in order of increasing p value. Don't know the nature of the estimates sent so safer to use p values
  #adding genes according to the order of estimates such as logFC or t-stat should not impact results much
  ordering = order(abs(genesCmap$pvalue), decreasing = FALSE)
  
  if(genesToStart >= 1)
  {
    while(genesToStart > 1)
      genesToStart = genesToStart/100
    
    warning("The genesToStart value supplied was greater than 1 and has been reduced to a fraction between 0 and 1")
    print(paste(genesToStart*100, "percent of the genes will be used in the first iteration of the analysis"))
    if(genesToStart == 1)
    {
      warning("As genesToSart was equal to 1 and would have crashed the program, it has been changed to 0.20")
      genesToStart = 0.20
    }
  }
  
  genesSigEnd = nrow(genesCmap)
  genesSigStart = ceiling(nrow(genesCmap)*genesToStart)
  genesSigInc = round((genesSigEnd - genesSigStart)/(numbIters - 1)) + 1
  
  print(paste0("Starting with ", genesSigStart, " genes and incrementing by ", genesSigInc, " genes. The last analysis will contain ", genesSigEnd," genes"))
  
  drugNotes = NULL
  if(!is.null(pharmSet))
    drugNotes = drugInfo(pharmSet)
  
  cmapList = list()
  pFrame = as.data.frame(NULL)
  print(paste("starting the",numbIters, "iterations"))
  for(g in 1:numbIters)
  {
    genesUse = genesSigStart + (g-1)*genesSigInc
    sigRows = ordering[1:genesUse]
    genesCmapUse = genesCmap[sigRows, ]
    
    cmapTab = as.data.frame(NULL)
    numDrugs = dim(drugPert)[2]
    drugNames = drugPert[1, ,1][1:numDrugs]
    drugNames = names(drugNames)
    noCores <- detectCores() - 2
    # Initiate cluster
    cl <- makeCluster(noCores)
    pco = library(PharmacoGx)
    #clusterEvalQ(cl, library(PharmacoGx))
    
    drugsList = list()
    for(i in 1:length(drugNameVec))
    {
      drugColInd = which(colnames(drugPert) == drugNameVec[i])
      if(drugEst == TRUE)
      {
        drugsList[[i]] = drugPert[,drugColInd ,c("tstat", "pvalue")]
      }else{
        drugsList[[i]] = drugPert[,drugColInd ,c("estimate", "pvalue")]
      }
    }
    
    csFunc = function(z) connectivityScore(x = genesCmapUse, y = z,method=drugScoreMeth, gwc.method=gwcMethod, nperm = numbPerms)
    environment(csFunc) <- .GlobalEnv
    
    clusterExport(cl, list("drugsList", "genesCmapUse", "numbPerms", "gwcMethod", "drugScoreMeth", "connectivityScore", "gwc", "intersectList", "pco", "corWeighted", "combineTest"), envir = environment())
    dataCor = clusterApplyLB(cl, drugsList, csFunc)
    stopCluster(cl)
    
    drugPVec = c()
    for(k in 1:length(drugNameVec))
    {
      drugPVec = c(drugPVec, dataCor[[k]]["p"])
    }
    pFrame = rbind(pFrame, drugPVec)
    rownames(pFrame)[nrow(pFrame)] = paste(genesUse, "genes")
    #write.table(cmapTabSort, paste("genes used is", genesUse, fileNameCmap), col.names = NA, sep="\t")
    print(g)
  }
  colnames(pFrame) = drugNameVec
  
  return(pFrame)
}

#library("devtools")
#library(roxygen2)
#setwd("C:\\Users\\micha\\Documents\\PMH Research\\GWC CMAP Package\\gwcCMap")
#document()
#setwd("..")
#install("gwcCmap")
#library(gwcCMap)

#setwd("C:\\Users\\micha\\Documents\\PMH Research\\Code from Neel Project")
#geneData = read.table("breast_normal_vs_dtp-all_cell_lines.txt",sep="\t", header=TRUE)
#geneIds = as.character(geneData$symbol)
#geneIds = as.character(geneData$ensembl_id)
#geneEsts = geneData$t
#pvals = geneData$P.Value
#volcPlotEsts = geneData$logFC

#setwd("C:\\Users\\micha\\Documents\\PMH Research\\GWC CMAP Package")
#load("cmap_sig_rna.RData")
#drugPert = drug.perturbation

#setwd("C:\\Users\\micha\\Documents\\PMH Research\\GWC CMAP Package")
#load("CMAP.RData")
#pharmSet = CMAP
#pCut=TRUE
#cutOff = .05
#genesToStart = 0.20
#numbIters = 10
#numbPerms = 1000
#gwcMethod = "spearman"
#drugEst = TRUE
#inspectResults = TRUE
#showMimic = TRUE
#mDataType = "rna"

#setwd("C:\\Users\\micha\\Documents\\PMH Research\\GWC CMAP Package")
#load("L1000_compounds.RData")
#pharmSet = L1000_compounds
#load("all_drugs.RData")
#drugPert = lincs_small.drugPerturbation

#drugList = rankDrugsGwc(geneIds, geneEsts, pvals, drugPert = drugPert[, 1:100, ], pharmSet = pharmSet, pCut = TRUE, cutOff = 0.05, genesToStart = 0.20, numbIters = 10, gwcMethod = "spearman", numbPerms = 1000, volcPlotEsts = volcPlotEsts, drugEst = TRUE, inspectResults = TRUE, showMimic = FALSE, mDataType = "rna")


