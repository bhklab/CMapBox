
#' Function to get drugs to use when benchmarking 
#'
#' This function will return a frame with drug clinical trial information for drugs that contained the search words inputted to this function in their study information. The function
#' is most useful for determining drugs to use for benchmarking in the benchMarkingCmap function. For instance, if one tests various
#' drug repurposing methods on a breast cancer signature, the input search word could be set to "breast" to find breast cancer drugs
#' to use when benchmarking
#' @param searchWords a character vector containing the key words to be used to select drug names
#' @param containsAll a boolean. If set to TRUE, a drug name is only returned if its clinical trial info contains all the search phrases/words. Default is FALSE
#' @param drugsToKeep a character vector containing the names of the drugs that should be searched by using the words supplied in searchWords. Used to prevent the returned drug names from containing drugs that wont be used in the benchmarking analysis as they are not in the drug perturbation signature used (ex send in CMap drug names)
#' @return a frame with clinical trial info for the drug names supplied, if available. If no info is available NA is returned
#' @export
#' @examples
#' searchWords = c("solid tumor", "depressive")
#' drugTrialFrame = getBenchmarkDrugs(searchWords)
#' #below line will show other studies the drugs have been involved in, studies that did not contain the search words
#' clinicalInfoFrame = getDrugClinicalInfo(drugTrialFrame$drug_name)

getBenchmarkDrugs = function(searchWords, containsAll = FALSE, drugsToKeep = NULL)
{
  data("cmapTrialInfo")
  
  if(!is.null(drugsToKeep))
  {
    drugsToKeep = tolower(drugsToKeep)
    repoNames = tolower(cmapTrialInfo$drug_name)
    
    overlapDrugInds = c()
    for(i in 1:length(drugsToKeep))
      overlapDrugInds = c(overlapDrugInds, which(repoNames == drugsToKeep[i]))
    cmapTrialInfo = cmapTrialInfo[overlapDrugInds, ] 
  }
  
  if(containsAll == FALSE){
    rowVec = c()
    for(i in 1:length(searchWords))
    {
      word = searchWords[i]
      for(j in 1:nrow(cmapTrialInfo))
        if(grepl(word, tolower(cmapTrialInfo$Indication[j])))
          rowVec = c(rowVec, j) 
    }
  }else{
    rowVec = c()
    for(i in 1:nrow(cmapTrialInfo))
    {
      addDrug = TRUE
      for(j in 1:length(searchWords))
      {
        word = searchWords[j]
        if(!grepl(word, tolower(cmapTrialInfo$Indication[i])))
          addDrug = FALSE
      }
      if(addDrug == TRUE)
        rowVec = c(rowVec, i) 
    }
  }
  drugInfo = cmapTrialInfo[unique(rowVec), ]
    
  return(drugInfo)
}
