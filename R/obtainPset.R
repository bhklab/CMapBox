
#' Function to download PSets for analysis (credit to PharmacoGx::downloadPset)
#'
#' The present package requires perturbation datasets (PharmacoSig objects with drug/gene interaction info) to obtain drug lists. These can be obtained
#' through this function which calls the pharmacoGx::downloadPSet function. Note that if the corresponding pharmacoSet object that was used to create
#' the pharmacoSig is also obtained it can be used to obtain additional plots that show the impact of the drugs on the individual genes via the seeDrugImpactGenes function
#' @param psetName The default is "CMAP" but "L1000_compounds" can also be used to obtain the L1000 PharmacoSig object
#' @return a PSet object
#' @export
#' @examples

obtainPset = function(psetName = "CMAP"){
  #update this function to obtain pharmSet object and drugPert objects. I.e get the pharmSet, use drugPert to get the pharmSig, and return both
  if(!(psetName == "CMAP" | psetName == "L1000_compounds"))
    stop(paste(psetName, "is not available, please supply either CMAP or L1000_compounds"))

  drugPert = PharmacoGx::downloadPSet(psetName)
  return(drugPert)
}

