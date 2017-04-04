
#' Function to get clinical trial info for drugs
#'
#' This function will return a data frame with clinical trial info for the drug names supplied, if available. Note that importing the data
#' file "cmapTrialInfo" will give one access to all of the available clinical trial info
#' @param drugNames a character vector containing the names of the drugs to return clinical trial info for
#' @return a frame with clinical trial info for the drug names supplied, if available. If no info is available NA is returned
#' @export
#' @examples
#' drugNames = c("doxorubicin", "cetuximab")
#' clinicalInfo = getDrugClinicalInfo(drugNames)

getDrugClinicalInfo = function(drugNames)
{
  data("cmapTrialInfo")
  rowInds = c()
  for(i in 1:length(drugNames))
  {
    rowInds = c(rowInds, which(tolower(cmapTrialInfo$drug_name) %in% tolower(drugNames[i])))
  }
  rowInds = which(tolower(cmapTrialInfo$drug_name) %in% tolower(drugNames))
  if(length(rowInds) > 0)
    clinicalInfo = cmapTrialInfo[rowInds, ]
  if(length(rowInds) == 0)
    clinicalInfo = NA
  return(clinicalInfo)
}
  

