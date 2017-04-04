#' Function to determine which genes to keep when there are multiples IDs per gene with different values
#'
#' This function returns the indices of genes in a vector that should be kept. Genes that should be kept are those without other genes having the same ID and genes with lower P values (default) or higher abs(estimates) of those genes with the same ID. 
#' @param geneIds a character vector containing the gene symbols, ensemble IDs, or entrez IDs for the genes to be analyzed
#' @param geneVals a numeric vector containing values or p values for the genes with the corresponding gene IDs
#' @param byP a boolean value specifying whether genes will be kept according to smallest P values or highest abs(value)
#' @return a numeric vector containing the indices for the geneId elements that should be kept
#' @export
#' @examples
#' data("geneDataGwc")
#' indsToKeep = getUniqueIds(geneIds = as.character(geneDataGwc$ensembl_id), geneVals = geneDataGwc$P.Value, byP = TRUE)
#' ##or keep genes with larger abs(logFC)
#' indsToKeep = getUniqueIds(geneIds = as.character(geneDataGwc$ensembl_id), geneVals = geneDataGwc$logFC, byP = FALSE)
#' geneDataUnique = geneDataGwc[indsToKeep, ]


getUniqueIds = function(geneIds, geneVals, byP)
{

  if(length(geneIds) != length(geneVals))
    warning("There is not a one to one correspondence between the gene IDs and gene values, an error will likely occur or the results will be incorrect")
  names(geneVals) = c(1:length(geneVals))
  names(geneIds) = c(1:length(geneIds))

  geneIdsPreSort = geneIds
  
  naInds = unique(c(which(is.na(geneIds))), c(which(is.na(geneVals))))
  if(length(naInds) > 0)
  {
    geneIds = geneIds[-naInds]
    geneVals = geneVals[-naInds]
  }
  sortInd = sort(geneIds, decreasing = FALSE, index.return=TRUE)$ix;
  #definitly sorts least to greatest
  
  geneIds = geneIds[sortInd];
  geneVals = geneVals[sortInd];
  #dataInfoProbeset = dataInfoProbeset[sortInd];
  #dataInfoGene = dataInfoGene[sortInd];
  
  bestProbes = c()
  i = 1;
  while(i < length(geneIds))
  {
    geneId = geneIds[i];
    origInd = i;
    probesWithId = c();
    while(geneIds[i] == geneId)
    {
      probesWithId = c(probesWithId, i);
      i = i + 1;
      if(i == length(geneVals))
      {
        if(geneIds[i] == geneId)
          probesWithId = c(probesWithId, i)
        if(geneIds[i] != geneId)
          bestProbes = c(bestProbes, i)
        
        break
      }
    }
    if(length(probesWithId) > 1)
    {
      #could add warning to let user know that if 2 probes that have the same ID have the same value the first is taken
      if(byP == TRUE)
      {
        probes = as.numeric(as.character(geneVals[probesWithId]))
        keepProbe = which(probes == min(probes))[1] + (origInd - 1);
      }else{
        probes = abs(as.numeric(as.character(geneVals[probesWithId])))
        keepProbe = which(probes == max(probes))[1] + (origInd - 1);
      }
      
      bestProbes = c(bestProbes, keepProbe)
    }
    if(length(probesWithId) == 1)
      bestProbes = c(bestProbes, origInd)
    
  }
  #get the indices of the probes to keep and then match them with the indices from the unsorted data
  keepRows = names(geneVals)[bestProbes]
  return(sort(as.numeric(keepRows)))
  
  
}


