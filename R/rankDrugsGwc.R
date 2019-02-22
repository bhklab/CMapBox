
#' Function to identify drug candidates.
#'
#' This function ranks a list of drugs according to their ability to reverse a phenotype. This is done by first computing scores for the drugs according to one of the various drug repurposing techniques over multiple trials, where each trial 
#' contains a different number of differentially expressed genes between the disease and normal phenotype. At each trial, the drugs are ranked relative to one another in terms of their ability to reverse the disease phenotype. Additionally, 
#' drugs are ranked according to how stable they are from trial to trial (as one varies the number of genes). A final score is given to each drug and is determined according to a drugs rank throughout each trial and its change in rank from trial to trial. For N drugs,
#' a final score of -N would imply the drug was the best at reversing the disease phenotype every single trial and that the drugs rank changed the least from one trial to the next, thus meaning it is both the best drug at reversing the disease phenotype and the most stable drug candidate.
#' @param geneIds a character vector containing the gene symbols, ensemble IDs, or entrez IDs for the genes to be analyzed
#' @param geneEsts a numeric vector containing estimates or directions (+1 -> up or -1 -> down) for the difference in expression between the two phenotyes. Note that positive values of this estimate (+1, or tstat, logFC, etc) must correspond to higher gene expression in the phenotype one would like to reverse
#' @param pvals P values from a t-test assessing the diferential expression of the genes between the two phenotypes. If not provided (i.e direction information only for geneEsts) then volcano plots for the data will not be generated and repeat gene IDs will be removed arbitrarily (first ID kept) as opposed to by p value
#' @param drugPert a drug perturbation signature (object of class PharmacoSig) with rownames that correspond to the ensemble IDs of the genes
#' @param pharmSet the PharmacoSet used to generate the drug perturbation signature. Supplying this will add additonal info about the drugs in the ranking tables.
#' @param drugScoreMeth a string specifying which drug repurposing technique to use to score the drugs during each trial. The options for this parameter are currently "gwc", "gwcCmapBox", "fgsea", and "xsum". Default is "gwcCmapBox"
#' @param pCut a boolean specifying whether to use pvalues to remove insignificant genes from the CMAP analysis (TRUE) or to remove genes from the CMAP analysis according to their supplied gene estimates (FALSE)
#' @param cutOff if pCut is TRUE then this value represents the p value threshold used to filter out genes. If pCut is FALSE then this value represents the fraction of genes present in both the data and drug perturbation signature with the top absolute value of gene estimates that will be left in the analysis. cutOff should be between 0 and 1. If p values ar enot provided, all genes supplied will be used in the analysis
#' @param genesToStart a value between 0 and 1 representing the fraction of genes present in the data and drug perturbation signature to use in the first iteration of the CMAP analysis. Recommended to be at least 0.40 to avoid large changes during the early iterations due to having too few genes present.
#' @param numbIters The number of iterations that will occur in the analysis. numbIters will set the rate at which genes are added to the analysis.
#' @param numbPerms The number of permutations to be used to compute the p value in the drug repurposing functions. 
#' @param gwcMethod a character string specifying which method to use when computing correlations in the gwc function. The options are spearman (default) or pearson.
#' @param volcPlotEsts a numeric vector containing gene estimates to be used for volcano plots, if not provided geneEsts will be used. If one desires to use t-stats for geneEsts in the gwc analysis it is recommended that one supply the logFC of the genes here if volcano plots are desired. 
#' @param drugEst a boolean specifying whether to use the estimates for each gene of the drug perturbation signature in the calculations (TRUE) or to use the t-stats for each gene in the drug perturbation signatur ein the calculations (FALSE). Default is TRUE
#' @param inspectResults a boolean specifying whether to display plots that will allow one to check that the correct drugs have been selected based on the data supplied. Default is TRUE (show plots)
#' @param showMimic a boolean (default is FALSE) specifying whether to show plots for the drug that upregulates and down regulates the genes that are overexpressed and under expressed in the disease state.
#' @param mDataType a string specifying the type of molecular data to retrieve molecular profiles for from the pharmacoSet, if one desires to inspect the results of the analysis (inspectResults = TRUE). Default is "rna"
#' @param extraData a data.frame where each column represents values one would like to inspect to determine if the genes corresponding to these values should be removed if they meet the conditions specified in extraCut and extraDirec. Useful if one would like to remove genes based on logFC or other values. extra data must have the same number of rows as the length of vectors geneIds, geneEsts, and pvals.
#' @param extraCut a numeric vector, with the nth value corresponding to the nth column of the m by n data.frame supplied in extraData, where each value indicates the value that needs to be reached for the columns of the data fram in order for the ids with that value to be kept or removed. Whether removal occurs when a value is greater or less than the value in this vector is specified in the extraDirec variable.
#' @param extraDirec a boolean vector, where TRUE means that ids whose value is greater than that specified in extraVals will be removed and FALSE means those with a value less than conditionVals will be removed. 
#' @return a data frame with information about the drugs in the analysis, including the drugs final scores. There are multiple connectivity scores, p values, and fdr adjusted p values in the results table as each value corresponds to one iteration in the analysis.
#' @keywords gwcCMap
#' @export
#' @examples
#' data("geneDataGwc")
#' data("drugPertEx")
#' data("psetSub")
#' drugResults = rankDrugsGwc(geneIds=geneDataGwc$ensembl_id, geneEsts=geneDataGwc$t, 
#'                            pvals = geneDataGwc$P.Value, drugPert = drugPertEx,
#'                            pharmSet=psetSub, drugScoreMeth = "gwc", 
#'                            volcPlotEsts = geneDataGwc$logFC)
#'                            
#' #below line shows how to additionally filter genes based on logFC
#' drugResults = rankDrugsGwc(geneIds=geneDataGwc$ensembl_id, geneEsts = geneDataGwc$t, 
#'                             pvals=geneDataGwc$P.Value, drugPert=drugPertEx, 
#'                             pharmSet=psetSub, volcPlotEsts=geneDataGwc$logFC, 
#'                             extraData = cbind(abs(geneDataGwc$logFC)), 
#'                             extraCut = c(0.5), extraDirec = c(TRUE))

rankDrugsGwc = function(geneIds, geneEsts, pvals = NULL, drugPert = NULL, pharmSet = NULL, drugScoreMeth = "gwcCmapBox", pCut = TRUE, cutOff = 0.05, genesToStart = 0.20, numbIters = 10, gwcMethod = "spearman", numbPerms = 1000, volcPlotEsts = NA, drugEst = TRUE, inspectResults = TRUE, showMimic = FALSE, mDataType = "rna", extraData = NA, extraCut = NA, extraDirec = NA)
{
  #########################################Section 1: Data Preperation#######################################################
  pvalsOrig = pvals
  if(is.null(pvals))
    pvals = sample(1:length(geneIds), length(geneIds))/length(geneIds)
  
  #convert geneEsts and pvals in case they are factors or character vectors
  geneEsts = as.numeric(as.character(geneEsts))
  pvals = as.numeric(as.character(pvals))
  
  geneDataDf = cleanData(geneIds, geneEsts, pvals, forRankAndPlot = TRUE)

  if(is.null(pvalsOrig) & nrow(geneDataDf) != length(geneIds))
    warning("repeat IDs were removed arbitrarily since p values were not supplied.")
  
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
    print("Computing perturbation signature from CMAP data.....this may take half a day, 
          considering finding and supplying a precomputed pset")
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
    if(pCut == TRUE & !is.null(pvalsOrig))
      print(paste("If it was intended that genes with a p value less than", cutOff, "should be kept, then the analysis will be conducted as intended"))
    #2 gene ests implies only direction is being used
    if(pCut == FALSE & length(unique(geneEsts)) > 2)
      print(paste("If it was intended that genes with an abs(estimate) in the top", cutOff*100, " percent should be kept, then the analysis will be conducted as intended"))
  }
  
  sigGenes = c(1:nrow(geneDataDf))
  if(pCut == TRUE & !is.null(pvalsOrig))
  {
    sigGenes = which(genesCmapAll$pvalue < cutOff)
    print(paste(length(sigGenes), "significant genes/genes with p <", cutOff,"in the analysis"))
  }else{
    if(length(unique(geneEsts)) > 2)
    {
      estOrd = order(abs(genesCmapAll$estimate), decreasing = TRUE)
      numKeep = floor(length(estOrd)*cutOff)
      sigGenes = estOrd[1:numKeep]
      print(paste(cutOff,"percent of genes with the top estimate will be used in the analysis for a total of", length(sigGenes), "genes"))
    }
  }
  
  genesCmap = genesCmapAll[sigGenes,]
  
  if(!is.na(extraData))
  {
    if(is.na(dim(extraData)))
    {
      extraData = cbind(extraData)
    }
    extraData = extraData[as.numeric(rownames(geneDataDf))]
    rowsToRem = removeGenes(idsOne = rownames(genesCmapAll), idsTwo = rownames(genesCmapAll), valuesTwo = extraData, conditionVals = extraCut, conditionsDirec = extraDirec)
    if(length(rowsToRem) > 0){
      genesCmap = genesCmap[-rowsToRem, ]
      print(paste("extraData was supplied and used to remove an additional ",length(rowsToRem) ," genes, the number of genes left in the analysis now is ", nrow(genesCmap)))
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
      warning("As genesToSart was equal to 1 and would have crashed the program, it has been changed to 0.40")
      genesToStart = 0.40
    }
  }
  
  genesSigEnd = nrow(genesCmap)
  genesSigStart = ceiling(nrow(genesCmap)*genesToStart)
  genesSigInc = round((genesSigEnd - genesSigStart)/(numbIters - 1)) + 1
  
  print(paste0("Starting with ", genesSigStart, " genes and incrementing by ", genesSigInc, " genes. The last analysis will contain ", genesSigEnd," genes"))
  
  drugNotes = NULL
  if(!is.null(pharmSet))
    drugNotes = drugInfo(pharmSet)
  
  #########################################Section 2: Drug Ranking###########################################################
  
  quickTest = FALSE
  cmapList = list()
  naDrugRowsVec = c()
  pFrame = as.data.frame(NULL)
  fdrFrame = as.data.frame(NULL)
  csFrame = as.data.frame(NULL)
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
    
    if(quickTest == TRUE)
    {
      numDrugs = 300
      numbPerms = 100 
    }
    
    drugsList = list()
    for(i in 1:numDrugs)
    {
      if(drugEst == TRUE)
      {
        drugsList[[i]] = drugPert[,i ,c("tstat", "pvalue")]
      }else{
        drugsList[[i]] = drugPert[,i ,c("estimate", "pvalue")]
      }
    }
  
    csFunc = function(z) connectivityScoreBroad(x = genesCmapUse, y = z,method=drugScoreMeth, gwc.method=gwcMethod, nperm = numbPerms)
    environment(csFunc) <- .GlobalEnv
    
    clusterExport(cl, list("drugsList", "genesCmapUse", "numbPerms", "drugScoreMeth", "gwcMethod", "connectivityScoreBroad", "gwc", "intersectList", "pco", "corWeighted", "combineTest"), envir = environment())
    dataCor = clusterApplyLB(cl, drugsList, csFunc)
    stopCluster(cl)
   
    for(i in 1:numDrugs)
    {
      if(!is.null(drugNotes))
      {
        drugInd = which(rownames(drugNotes) == drugNames[i])
        drugInfoTab = drugNotes[drugInd, ]
        cmapTab = rbind(cmapTab, c(dataCor[[i]], as.character(drugInfoTab)), stringsAsFactors = FALSE) 
        colnames(cmapTab) = c("score", "pvalue", colnames(drugNotes[1,]))
      }else{
        cmapTab = rbind(cmapTab, c(dataCor[[i]]), stringsAsFactors = FALSE) 
        colnames(cmapTab) = c("score", "pvalue")
      }
    }
    
    rownames(cmapTab) = drugNames[1:numDrugs]
    #sorting = order(cmapTab[,"score"], decreasing = TRUE)
    #cmapTabSort = cmapTab[sorting,]
    cmapTabNew = cbind(cmapTab[, 1:2], p.adjust(cmapTab[,2], method = "fdr"), cmapTab[, 3:ncol(cmapTab)])
    colnames(cmapTabNew)[3] = "fdr adjusted p"
    cmapTab = cmapTabNew

    
    #one of the drugs score kept showing up as NA, so remove it
    naDrugRows = which(is.na(cmapTab[,1]))
    naDrugRowsVec = c(naDrugRowsVec, naDrugRows)
    #if(length(naDrugRows) > 0)
    #  cmapTab = cmapTab[-naDrugRows, ]
    
    #rbind(paste(as.character(vec), collapse = ", "))
    if(g == 1)
    {
      pFrame = cbind(format(as.numeric(as.character(cmapTab[, "pvalue"])), scientific =  TRUE, digits = 3))
      fdrFrame = cbind(format(as.numeric(as.character(cmapTab[, "fdr adjusted p"])), scientific =  TRUE, digits = 3))
      csFrame = cbind(sprintf("%.5f", as.numeric(as.character(cmapTab[, "score"]))))
    }else{
      pFrame = cbind(pFrame, format(as.numeric(as.character(cmapTab[, "pvalue"])), scientific =  TRUE, digits = 3))
      fdrFrame = cbind(fdrFrame, format(as.numeric(as.character(cmapTab[, "fdr adjusted p"])), scientific =  TRUE, digits = 3))
      csFrame = cbind(csFrame, sprintf("%.5f", as.numeric(as.character(cmapTab[, "score"]))))
    }
      
    
    cmapList[[g]] = cmapTab
    #write.table(cmapTabSort, paste("genes used is", genesUse, fileNameCmap), col.names = NA, sep="\t")
    print(g)
  }
  rownames(pFrame) = drugNames[1:numDrugs]
  rownames(fdrFrame) = drugNames[1:numDrugs]
  rownames(csFrame) = drugNames[1:numDrugs]
  
  if(length(naDrugRowsVec) > 0)
  {
    naDrugRowsVec = unique(naDrugRowsVec)
    pFrame = pFrame[-naDrugRowsVec, ]
    csFrame = csFrame[-naDrugRowsVec, ]
    fdrFrame = fdrFrame[-naDrugRowsVec, ]
    for(g in 1:numbIters)
    {
      cmapList[[g]] = cmapList[[g]][-naDrugRowsVec, ]
    }
  }
  #rank the drugs according to their connectivity score and compute the change in their rank over the iterations
  rankFrame = as.data.frame(NULL)
  #below frame finds drugs that mimick the phenotype the most
  rankFrameRev = as.data.frame(NULL)
  rankChangeFrame = as.data.frame(NULL)
  colnameVec = c()
  for(g in 1:numbIters)
  {
    ranking = c()
    rankingRev = c()
    cmapRes = cmapList[[g]]["score"]
    ordering = sort(-1*as.numeric(cmapRes[,1]), index.return=TRUE,decreasing = TRUE)$ix
    orderingRev = sort(as.numeric(cmapRes[,1]), index.return=TRUE,decreasing = TRUE)$ix
    for(i in 1:length(ordering))
    {
      ranking[ordering[i]] = i
      rankingRev[orderingRev[i]] = i
    }
    if(g == 1)
    {
      rankFrame = as.data.frame(ranking) 
      rankFrameRev = as.data.frame(rankingRev)
    }
    if(g > 1)
    {
      rankFrame = cbind(rankFrame, ranking) 
      rankFrameRev = cbind(rankFrameRev, rankingRev)
    }
    if(g == 2)
      rankChangeFrame = as.data.frame(abs(ranking - oldRanking))
    if(g > 2)
      rankChangeFrame = cbind(rankChangeFrame, abs(ranking - oldRanking))
    colnameVec = c(colnameVec, as.character(genesSigStart + (g-1)*genesSigInc))
    oldRanking = ranking
  }
  rownames(rankFrame) = rownames(cmapRes)
  colnames(rankFrame) = colnameVec
  
  rownames(rankFrameRev) = rownames(cmapRes)
  colnames(rankFrameRev) = colnameVec
  
  rownames(rankChangeFrame) = rownames(cmapRes)
  colnames(rankChangeFrame) = colnameVec[2:length(colnameVec)]
  
  rankChangeFrameOrig = rankChangeFrame
  for(g in 1:dim(rankChangeFrame)[2])
  {
    changes = rankChangeFrame[,g]
    changeRanks = rank(-1*changes, ties.method = "min")
    changeRanksOrig = rank(changes, ties.method = "min")
    rankChangeFrame[,g] = changeRanks
    rankChangeFrameOrig[,g] = changeRanksOrig
  }
  
  #if a drug ranked higher on the list that ranked drugs according to negative connectivity scores use the score from that list as its final
  #score and multiply the score by a negative 1. Otherwise use the rank from other (positive connectivity score) list
  drugRankMeansPos = rowMeans(rankFrame)
  drugRankMeansNeg = rowMeans(rankFrameRev)
  drugRankSigns = c()
  drugRankMeans = c()
  for(g in 1:length(drugRankMeansPos))
  {
    if(drugRankMeansPos[g] >= drugRankMeansNeg[g])
    {
      drugRankMeans[g] = drugRankMeansPos[g]
      drugRankSigns[g] = 1
    }
    if(drugRankMeansNeg[g] > drugRankMeansPos[g])
    {
      drugRankMeans[g] = drugRankMeansNeg[g]
      drugRankSigns[g] = -1
    }
  }

  drugRankChangeMeans = rowMeans(rankChangeFrame)
  drugRankChangeMeansOrig = rowMeans(rankChangeFrameOrig)
  #final score is an average of the drugs mean rank and the rank of its stability from iteration to iteration
  sigFinalScores = drugRankSigns*(drugRankMeans + drugRankChangeMeans)/2
  
  pCol = c()
  fdrCol = c()
  csCol = c()
  for(i in 1:nrow(csFrame))
  {
    pCol[i] = paste(pFrame[i, ], collapse = ", ")
    fdrCol[i] = paste(fdrFrame[i, ], collapse = ", ")
    csCol[i] = paste(csFrame[i, ], collapse = ", ")
  }
  #although best final scores for reversing the phenotypes are large negative ranks, paste drugRankMeansPos and drugRankMeansOrig to the table
  #as they are easiest to interpret. A value of 1 on these lists means the drug had the largest negative connectivity score each iteration (rankMeans)
  #and its rank changed the least from iteration to iteration relative to all the other drugs (rankChange)
  cmapResults = cbind(sigFinalScores, drugRankMeansPos, drugRankChangeMeansOrig, csCol, pCol, fdrCol, cmapList[[1]][4:ncol(cmapList[[1]])])
  sorting = order(cmapResults[,"sigFinalScores"], decreasing = FALSE)
  cmapResults = cmapResults[sorting,]
  colnames(cmapResults)[1:6] = c("Final Score", "Mean Rank (1 implies most negative connectivity scores)", "Rank of Mean Rank Change (lower implies more stable)", "connectivity scores", "p values", "fdr adjusted p")
  
  #rankFrameInd = which(rownames(rankFrame) == topDrug)
  #rankChangeInd = which(rownames(rankChangeFrame) == topDrug)
  #rankVsSize = rankFrame[rankFrameInd, ]
  #changeVsSize = rankChangeFrameOrig[rankChangeInd,]
  
  #if(sum(grepl("BRD-", rownames(cmapResults))) > nrow(cmapResults)/2)
  #{
  #  rownames(cmapResults) = paste(cmapResults$pert_iname, " (ID ", rownames(cmapResults),")", sep = "")
  #  #give drugPert appropr
  #  for(i in 1:length(colnames(drugPert)))
  #    colnames(drugPert)[i] = paste(drugNotes$pert_iname[which(colnames(drugPert)[i] == rownames(drugNotes))], " (ID ",colnames(drugPert)[i],")", sep = "")
      
  #}
  #topDrug = rownames(cmapResults)[1]
  
#########################################Section 3: Results Validation###########################################################
  
  #gsea ranks vs genes seem to trend upwards, add genes previous best ranked drug has higher ranks, gwc seems to stabilize
  if(inspectResults == TRUE)
  {
    #sizeVec = c()
    #for(g in 1:numbIters)
    #  sizeVec = c(sizeVec, genesSigStart + (g-1)*genesSigInc)
    
    #plot(sizeVec, rankVsSize, main = paste("Rank as a Function of Query Size for the Top Drug", topDrug), xlab = "Query Size (number of genes in the analysis)", ylab = "Rank (1 implies best at reversing disease phenotype)")
    #plot(1:(numbIters-1), changeVsSize, main = paste("Change in Rank from Query to Query for the Top Drug", topDrug), xlab = "Query Number", ylab = "abs(Rank of Mean Rank Change) Between Queries")
    
    direcOnly = FALSE
    if(is.null(pvalsOrig) | length(unique(geneEsts)) == 2)
      direcOnly = TRUE
    
    topDrug = rownames(cmapResults)[1]
    drugRankPlots(cmapResults, topDrug)
    
    if(is.na(volcPlotEsts[1]) | is.null(pvalsOrig)){
      volcPlotEsts = geneDataDf$geneEsts
    }else{
      volcPlotEsts = volcPlotEsts[as.numeric(rownames(geneDataDf))]
    }

    idsUsed = rownames(genesCmap)
    idsUsed = which(geneDataDf$ensemble %in% idsUsed)
    inspectDrugResults(volcPlotEsts[idsUsed], geneDataDf$symbol[idsUsed], geneDataDf$pvals[idsUsed], geneDataDf$ensemble[idsUsed], drugPert, pharmSet, drugScoreMeth, mDataType = mDataType, cmapResults, gwcMethod, drugEst, direcOnly, drugVolcPlot = FALSE, showMimic = showMimic)
  }
  
  #add clinical trial info if available
  data("cmapTrialInfo")
  cmapNames = tolower(rownames(cmapResults))
  repoNames = tolower(cmapTrialInfo$drug_name)

  if(length(which(cmapNames %in% repoNames)) > 0)
  {
    trialInfoVec = c()
    trialStatusVec = c()
    for(i in 1:nrow(cmapResults))
    {
      trialInds = which(tolower(rownames(cmapResults)[i]) == tolower(cmapTrialInfo$drug_name))
      if(length(trialInds) > 0){
        drugClinInfo = cmapTrialInfo[trialInds, ]
        trialInfo = drugClinInfo$Indication
        trialInfoVec[i] = paste(trialInfo, collapse = " | ")
        trialStatus = drugClinInfo$status
        trialStatusVec[i] = paste(trialStatus, collapse = " | ")
      }else{
        trialInfoVec[i] = NA
        trialStatusVec[i] = NA
      }
    } 
    cmapResults = cbind(cmapResults, trialInfoVec, trialStatusVec)
    colnames(cmapResults)[(ncol(cmapResults) - 1):(ncol(cmapResults))] = c("Drug Clinical Trials", "Trial Outcomes")
  }
  #double check that each method returns connectivity scores between 0 and 1!
  
  return(cmapResults)
}
#to add new drug repurposing technique
#1. add technique to connectivityScore function and ensure more negative outcome implies reverses the phenotype (for ranking later)
#2. adjust getSigGeneVolcPlot to be able to determine which genes drove the connectivity score using that technique
#3. add statement to inspectDrugResults function explaining how to interpret data to determine if the analysis is correct
#4. adjust drugScoreMeth in file info section to include the technique
#5. adjust vignette or other areas that should make reference to the technique

#1250 genes at iteration 10 in original data received

#library("devtools")
#library(roxygen2)
#setwd("C:\\Users\\micha\\Documents\\PMH Research\\GWC CMAP Package\\gwcCMap")
#document()
#setwd("..")
#install("gwcCmap")
#library(gwcCMap)
