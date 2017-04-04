#' Function to remove genes from a dataset based on specified conditions
#'
#' This function removes genes from data set one from data set two according to the supplied conditions on the values sent to the function
#' @param idsOne a character vector containing the gene ids from the data that one would like to remove ids from based on data set two 
#' @param idsTwo a character vector containing the genes ids for the data supplied that one would like to use to determine which ids to remove from data set one
#' @param valuesTwo a data.frame where each column represents values one would like to inspect to determine if the ids corresponding to these values should be removed if they meet a specified condition
#' @param conditionVals a numeric vector, with the nth value corresponding to the nth column of the m by n data.frame supplied, where each value indicates the value that needs to be reached for the columns of the data fram in order for the ids with that value to be removed. Whether removal occurs when a value is greater or less than the value in this vector is specified in the conditionsDirec variable.
#' @param conditionsDirec a boolean vector, where TRUE means that ids whose value is greater than that specified in conditionVals will be removed and FALSE means those with a value less than conditionVals will be removed. 
#' @return the indices of ids in idsOne that should be removed as they dont meet the specified criteria
#' @export
#' @examples
#' data("geneDataGwc")
#' geneDataClean = cleanData(geneIds = geneDataGwc$ensembl_id, geneEsts = geneDataGwc$logFC, pvals = geneDataGwc$P.Value)
#' #keep genes with p value < 0.05 and |logFC| > .5
#' rowsToRem = removeGenes(idsOne = geneDataClean$ensemble, idsTwo = geneDataClean$ensemble, valuesTwo = cbind(geneDataClean$pvals, abs(geneDataClean$geneEsts)), conditionVals = c(0.05, .5), conditionsDirec = c(TRUE, FALSE))
#' sigGeneData = geneDataClean[-rowsToRem, ]
#' 


removeGenes = function(idsOne, idsTwo, valuesTwo, conditionVals, conditionsDirec)
{
  if(is.null(dim(valuesTwo)))
  {
    valuesTwo = cbind(valuesTwo)
  }
  rmInds = rep(TRUE, length(idsTwo))
  #rmLast = rep(length(idsTwo), TRUE)
  for(i in 1:ncol(valuesTwo))
  {
    if(conditionsDirec[i] == TRUE)
    {
      rmCur = valuesTwo[, i] > conditionVals[i]
    }else{
      rmCur = valuesTwo[, i] < conditionVals[i]
    }
    rmInds = rmCur & rmInds
    #rmInds = rmCur & rmLast
    rmLast = rmInds
  }
  rmInds = which(rmInds == TRUE)
  idsRm = idsTwo[rmInds]
  idIndsInDataOne = which(idsOne %in% idsRm)
  
  #below is equivalent, above just returns indices out of order
  #idIndsInDataOne = c()
  #for(i in 1:length(idsRm))
  #  idIndsInDataOne = c(idIndsInDataOne, which(idsOne == idsRm[i]))
  
  return(idIndsInDataOne)
}


