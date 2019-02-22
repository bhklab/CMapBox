#' Function to compute drug ranks and plot their rank and change in rank over the iterations/as the number of genes in the analysis change 
#'
#' This function uses the results from rankDrugsGwc to plot the rank of a drug over the iterations and the rank of the change in rank of the drug over the iterations (to see drug stability over the iterations)
#' @param drugResults a data frame returned from rankDrugsGwc containing the connectivity scores of the drugs over the iterations
#' @param drugName The name of the drug one would like to generate plots for to examine its performance as a function of the number of genes in the analysis. Defaults to the drug with the most negative (most likely to reverse disease phenotype) rank if a name is not provided.
#' @keywords gwcCMap
#' @export
#' @examples
#' data("geneDataGwc")
#' data("drugPertEx")
#' geneDataClean = cleanData(geneIds = geneDataGwc$symbol, geneEsts = geneDataGwc$t, pvals = geneDataGwc$P.Value)
#' pFrame = getDrugsPValue(drugNameVec = c("BRD-K78431006", "BRD-A68739437", "BRD-K73368362"), geneIds = geneDataClean$ensemble, geneEsts = geneDataClean$geneEsts, pvals = geneDataClean$pvals, drugPert = drugPertEx, numbIters = 5, numbPerms = 50000)

drugRankPlots = function(drugResults, drugName = NULL)
{
  if(is.null(drugName))
    drugName = rownames(drugResults)[1]
  
  #recreate cmapList object from rankDrugsGwc and then reuse the code from that file
  cmapList = list()
  numbIters = length(as.numeric(strsplit(as.vector(drugResults$`connectivity scores`[1]), ",")[[1]]))
  for(i in 1:numbIters)
  {
    cmapTab = as.data.frame(NULL)
    scores = c()
    for(j in 1:nrow(drugResults))
    {
      scores = c(scores, as.numeric(strsplit(as.vector(drugResults$`connectivity scores`[j]), ",")[[1]])[i])
    }
    cmapTab = cbind(scores)
    rownames(cmapTab) = rownames(drugResults)
    colnames(cmapTab) = c("score")
    cmapList[[i]] = as.data.frame(cmapTab)
  }
  
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
    #colnameVec = c(colnameVec, as.character(genesSigStart + (g-1)*genesSigInc))
    oldRanking = ranking
  }
  rownames(rankFrame) = rownames(cmapRes)
  #colnames(rankFrame) = colnameVec
  
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
    if(drugRankMeansPos[g] > drugRankMeansNeg[g])
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
  for(i in 1:nrow(pFrame))
  {
    pCol[i] = paste(pFrame[i, ], collapse = ", ")
    fdrCol[i] = paste(fdrFrame[i, ], collapse = ", ")
    csCol[i] = paste(csFrame[i, ], collapse = ", ")
  }
  
  #although best final scores for reversing the phenotypes are large negative ranks, 
  #paste drugRankMeansPos and drugRankMeansOrig to the table
  #as they are easiest to interpret. A value of 1 on these lists means the drug 
  #had the largest negative connectivity score each iteration (rankMeans)
  #and its rank changed the least from iteration to iteration relative to all 
  #the other drugs (rankChange)
  cmapResults = cbind(sigFinalScores, drugRankMeansPos, drugRankChangeMeansOrig, csCol, pCol, fdrCol)
  #sorting = order(cmapResults[,"sigFinalScores"], decreasing = FALSE)
  #cmapResults = cmapResults[sorting,]
  colnames(cmapResults)[1:6] = c("Final Score", "Mean Rank (1 implies most negative connectivity scores)", "Rank of Mean Rank Change (lower implies more stable)", "connectivity scores", "p values", "fdr adjusted p")
  
  topDrug = drugName
  rankFrameInd = which(rownames(rankFrame) == topDrug)
  rankChangeInd = which(rownames(rankChangeFrame) == topDrug)
  rankVsSize = rankFrame[rankFrameInd, ]
  changeVsSize = rankChangeFrameOrig[rankChangeInd,]
  
  plot(1:length(rankVsSize), rankVsSize, main = paste("Rank as a Function of the Iterations for", topDrug), xlab = "Query Number", ylab = "Rank (1 implies best at reversing disease phenotype)")
  plot(1:length(changeVsSize), changeVsSize, main = paste("Change in Rank from Query to Query for", topDrug), xlab = "Query Number", ylab = "abs(Rank of Mean Rank Change) Between Queries")
}

