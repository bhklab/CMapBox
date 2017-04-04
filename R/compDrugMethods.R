
#' Function to compare drug repurposing techniques.
#'
#' This function operates like the rankDrugsGwc function but allows users to compare multiple drug repurposing methods. Users can either compare different methods,
#' the results of using the same method with different parameters, or can run the packges main function on different data sets. This is accomplished by passing in 
#' lists instead of vectors or values, where the n'th list element is used in the n'th analysis. Note that if one is interested in keeping the variables the same but trying different methods,
#' only the drugScoreMeth variable would be in list format. If some variables are lists and other are not or are shorter lists, the supplied variable or last list element
#' will be used in the analyses containing the longer lists variables. A data frame is returned that contains the ranks of the drugs in each method 
#' and the average rank of the drug across mutliple methods.
#' @param numAnalyses an integer specifying how many analyses are intended/the max size of one of the lists passed to the function. Note that not providing 
#' @param geneIds a list of character vectors containing the gene symbols, ensemble IDs, or entrez IDs for the genes to be analyzed in each analysis
#' @param geneEsts a list of numeric vectors containing estimates for the difference in expression between the two phenotyes in each analysis. Note that positive values of this estimate (tstat, logFC, etc) must correspond to higher gene expression in the phenotype one would like to reverse
#' @param pvals a list of numeric vectors containing P values from a t-test assessing the diferential expression of the genes between the two phenotypes in each analysis
#' @param drugPert a drug perturbation signature (object of class PharmacoSig) with rownames that correspond to the ensemble IDs of the genes
#' @param pharmSet the PharmacoSet used to generate the drug perturbation signature. Supplying this will add additonal info about the drugs in the ranking tables.
#' @param drugScoreMeth a list of strings specifying which drug repurposing technique to use to score the drugs during each trial in an analysis. The options for this parameter are currently "gwc" and "fgsea". Default is "gwc"
#' @param pCut a list of booleans for each analysis specifying whether to use pvalues to remove insignificant genes from the analysis (TRUE) or to remove genes from the analysis according to their supplied gene estimates (FALSE)
#' @param cutOff if pCut is TRUE then this value represents the p value threshold used to filter out genes. If pCut is FALSE then this value represents the fraction of genes present in both the data and drug perturbation signature with the top absolute value of gene estimates that will be left in the analysis. cutOff should be between 0 and 1.
#' @param genesToStart a list of values between 0 and 1 representing the fraction of significant genes present in the data and drug perturbation signature to use in the first iteration of each analysis. Recommended to be at least 0.15 to avoid large changes during the early iterations due to having too few genes present.
#' @param numbIters a list of integers specifying the number of iterations that will occur in each analysis. numbIters will set the rate at which genes are added to the analysis.
#' @param numbPerms a list of integers specifying the number of permutations to be used to compute the p value in the drug repurposing methods which use permutation tests to calculate p values. 
#' @param gwcMethod a list of strings specifying which method to use when computing correlations in the gwc function in each analysis (if gwc is being used). The options are spearman (default) or pearson.
#' @param volcPlotEsts a list of numeric vectors containing gene estimates to be used for volcano plots, if not provided geneEsts will be used. However, if one desires to use t-stats for geneEsts in the gwc analysis it is recommended that one use another estimate (eg logFC) if volcano plots are desired due to the high correlation between t-stats and p values. 
#' @param drugEst a list of booleans for each analysis specifying whether to use the estimates for each gene of the drug perturbation signature in the gwc calculation (TRUE) or to use the t-stats for each gene in the drug perturbation signature in the calculations (FALSE). Default is TRUE
#' @param extraData a list of data.frames to be use din each analysis, where each column represents values one would like to inspect to determine if the genes corresponding to these values should be removed if they meet the conditions specified in extraCut and extraDirec. Useful if one would like to remove genes based on logFC or other values. extra data must have the same number of rows as the length of vectors geneIds, geneEsts, and pvals.
#' @param extraCut a list of numeric vectors, with the nth value corresponding to the nth column of the m by n data.frame supplied in extraData, where each value indicates the value that needs to be reached for the columns of the data fram in order for the ids with that value to be kept or removed. Whether removal occurs when a value is greater or less than the value in this vector is specified in the extraDirec variable.
#' @param extraDirec a list of boolean vectors, where TRUE means that ids whose value is greater than that specified in extraVals will be removed and FALSE means those with a value less than conditionVals will be removed. 
#' @return a data frame with information about the drugs in each analysis, including the drugs final scores and the drugs average final scores across the different methods. 
#' @keywords gwcCMap
#' @export
#' @examples
#' data("geneDataGwc")
#' data("drugPertEx")
#' data("psetSub")
#' #compare a gwc, fgsea, and xsum, analysis that have the same parameters and also a gwc analysis with different parameters
#' compareResults = compDrugMethods(4, geneIds = geneDataGwc$ensembl_id, geneEsts = geneDataGwc$t, pvals = geneDataGwc$P.Value, drugPert = drugPertEx, pharmSet = psetSub, drugScoreMeth = list("gwc", "fgsea","xsum", "gwc"), pCut = list(TRUE, TRUE, TRUE, FALSE), cutOff = list(0.05, 0.05, 0.05, 0.10), genesToStart = list(0.20, 0.20, 0.20, 0.10), numbIters = 5, gwcMethod = "spearman", numbPerms = 1000, volcPlotEsts = NULL, drugEst = TRUE, extraData = NA, extraCut = NA, extraDirec = NA)


compDrugMethods = function(numAnalyses, geneIds, geneEsts, pvals, drugPert, pharmSet, drugScoreMeth, pCut = TRUE, cutOff = 0.05, genesToStart = 0.20, numbIters = 10, gwcMethod = "spearman", numbPerms = 1000, volcPlotEsts = NULL, drugEst = TRUE, extraData = NA, extraCut = NA, extraDirec = NA)
{
  #create list of lists and pass params to rankDrugsGwc
  #could use match.call to make it cleaner but would be less user friend
  paramListsOrig = list(numAnalyses, geneIds, geneEsts, pvals, drugPert, pharmSet, drugScoreMeth, pCut, cutOff, genesToStart, numbIters, gwcMethod, numbPerms, volcPlotEsts, drugEst, extraData, extraCut, extraDirec)
  numbArgs = length(paramListsOrig)
  paramList = list()
  tempList = list()
  for(i in 1:numbArgs)
  {
    if(class(paramListsOrig[[i]]) == "list"){
      for(j in 1:length(paramListsOrig[[i]]))
      {
        tempList[[j]] = paramListsOrig[[i]][[j]]
      }
      if(length(tempList) < numAnalyses)
      {
        lastListInd = length(tempList)
        for(j in (length(tempList+1):numAnalyses))
          tempList[[j]] = tempList[[lastListInd]]
      }
    }else if(class(paramListsOrig[[i]]) == "NULL"){
      for(j in 1:numAnalyses)
      {
        tempList[[j]] = NA
      }
    }else{
      for(j in 1:numAnalyses)
      {
        tempList[[j]] = paramListsOrig[[i]]
      }
    }
    paramList[[i]] = tempList
  }
  #length should be one for paramList
  #7,8,9,10,11,12,13,15,16,17,18
  listElsOne = c(7,8,9,10,11,12,13,15,16,17,18)
  for(i in 1:length(listElsOne))
  {
    if(length(paramList[[listElsOne[i]]][1]) > 1)
    {
      stop("it appears that you supplied a vector for drugScorMeth (ex c(\'\'gwc\'\', \'\'fgsea\'\') instead of a list(\'\'gwc\'\', \'\'fgsea\'\'). 
           please rerun the function passing lists in to separate vectors and variables that should differ between each analysis and not vectors.")
    } 
  }
  drugResList = list()
  for(i in 1:numAnalyses)
  {
    drugResList[[i]] = rankDrugsGwc(paramList[[2]][[i]], paramList[[3]][[i]], paramList[[4]][[i]], paramList[[5]][[i]], paramList[[6]][[i]], paramList[[7]][[i]], paramList[[8]][[i]], paramList[[9]][[i]], paramList[[10]][[i]], paramList[[11]][[i]], paramList[[12]][[i]], paramList[[13]][[i]], paramList[[14]][[i]], paramList[[15]][[i]], inspectResults = FALSE, showMimic = FALSE, mDataType = "rna", paramList[[16]][[i]], paramList[[17]][[i]], paramList[[18]][[i]])
    drugResList[[i]] = drugResList[[i]][order(rownames(drugResList[[i]])), ]
  }
  
  compRankFrame = cbind(rep(0, nrow(drugResList[[1]])))
  for(i in 1:numAnalyses)
  {
    compRankFrame = cbind(compRankFrame, drugResList[[i]]$`Final Score`)
    colnames(compRankFrame)[ncol(compRankFrame)] = paste("Final Score", drugScoreMeth[[i]], "- Analysis", i)
  }
  compRankFrame[, 1] = rowSums(compRankFrame[, 2:ncol(compRankFrame)])/(ncol(compRankFrame) - 1)
  colnames(compRankFrame)[1] = "Final Score Averaged Over the Analyses"
  compRankFrame = cbind(compRankFrame, drugResList[[1]][, 7:ncol(drugResList[[1]])])
  sorting = order(compRankFrame[,1], decreasing = FALSE)
  compRankFrame = compRankFrame[sorting,]
  return(compRankFrame)
}

# geneIds = geneDataGwc$ensembl_id
# geneEsts = geneDataGwc$t
# pvals = geneDataGwc$P.Value
# drugPert = drugPertEx
# pharmSet = psetSub
# drugScoreMeth = list("gwc", "fgsea")
# volcPlotEsts = geneDataGwc$logFC
# numAnalyses = 1
# pCut = TRUE
# cutOff = 0.05
# genesToStart = 0.20
# numbIters = 5
# gwcMethod = "spearman"
# numbPerms = 1000
# volcPlotEsts = NA
# drugEst = TRUE
# extraData = NA
# extraCut = NA
# extraDirec = NA

