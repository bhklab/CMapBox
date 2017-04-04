#' Function to compute a gene set enrichment analysis on a drug in a pharmacoSet
#'
#' This function removes genes from data set one from data set two according to the supplied conditions on the values sent to the function
#' @param pathways list of gene sets to check. Gene IDs in the gene sets must be entrez IDs or ensemble IDs
#' @param drugName the name of the drug in the pharmacoSet to perform the analysis on
#' @param drugPert a drug perturbation signature (object of class PharmacoSig) with rownames that correspond to either entrez IDs or ensemble IDs for the genes.
#' @param useTstat a boolean specifying whether to use the drugs t-stat in fgsea (default is TRUE). If FALSE the drug's estimate for each gene will be used
#' @param nperm number of permutations to do, default is 50000
#' @return A table with GSEA results. Each row corresponds to a tested pathway. See fgsea for info pertaining to the columns
#' @export
#' @examples
#' data("drugPertEx")
#' drugGseaResults = runFgseaPharmDrug(pathways = "cc", drugName = "BRD-K78431006", drugPert = drugPertEx)

runFgseaPharmDrug = function(pathways, drugName, drugPert, useTstat = TRUE, nperm = 50000)
{
  entIds = FALSE
  
  if(is.element("GSA", installed.packages()[,1]) == FALSE)
    install.packages("GSA")
  library(GSA)
  if(is.element("piano", installed.packages()[,1]) == FALSE)
    install.packages("piano")
  library(piano)
  if(is.element("parallel", installed.packages()[,1]) == FALSE)
    install.packages("parallel")
  library(parallel)
  if(is.element("fgsea", installed.packages()[,1]) == FALSE)
    install.packages("fgsea")
  library(fgsea)
  
  numCores = 1
  if((detectCores() - 1) > 0)
    coresNum = detectCores() - 1
  
  pathwaysObj = pathways
  if(pathwaysObj[1] == "all")
  {
    entIds = TRUE
    data("aGwcPack")
    pathwaysObj = aGwcPack
  }else if(pathwaysObj[1] == "bp"){
    entIds = TRUE
    data("abpGwcPack")
    pathwaysObj = abpGwcPack
  }else if(pathwaysObj[1] == "cc"){
    entIds = TRUE
    data("accGwcPack")
    pathwaysObj = accGwcPack
  }else if(pathwaysObj[1] == "mf"){
    entIds = TRUE
    data("amfGwcPack")
    pathwaysObj = amfGwcPack
  }
  
  if(dim(drugPert)[3] > dim(drugPert)[2])
    drugPert = aperm(drugPert, c(1,3,2))
  #deal with CMAP rownames _at
  if(grepl("_at", rownames(drugPert)[1]))
  {
    rownames(drugPert) = gsub("_at_", "", paste0(rownames(drugPert), "_"))
    entIds = FALSE
  }
  #deal with L1000 rownames geneid.
  if(grepl("geneid.", rownames(drugPert)[1]) == TRUE)
  {
    geneEntrezs = gsub("geneid.", "", rownames(drugPert))
    mapping = AnnotationDbi::select(org.Hs.eg.db, keytype="ENTREZID", keys=as.character(geneEntrezs), columns=c("ENSEMBL"))
    mapping = mapping[!duplicated(mapping[,1]),]
    rownames(drugPert) = mapping[, 2]
  }
  if(grepl("ENSG", rownames(drugPert)[1]) == FALSE)
  {
    geneEntrezs = rownames(drugPert)
    mapping = AnnotationDbi::select(org.Hs.eg.db, keytype="ENTREZID", keys=as.character(geneEntrezs), columns=c("ENSEMBL"))
    mapping = mapping[!duplicated(mapping[,1]),]
    rownames(drugPert) = mapping[, 2]
  }
  #now drugPert/drugGeneStats rownames/names are in ensemble ids
  
  drugPertInd = which(colnames(drugPert) == drugName)
  if(useTstat == TRUE)
  {
    drugGeneStats = drugPert[,drugPertInd ,c("tstat")] 
  }else{
    drugGeneStats = drugPert[,drugPertInd ,c("estimate")]
  }
  
  #entIds true implies user has asked for one of the preset pathway objects
  if(entIds == TRUE){
    mapping = AnnotationDbi::select(org.Hs.eg.db, keytype="ENSEMBL", keys=as.character(names(drugGeneStats)), columns=c("ENTREZID"))
    mapping = mapping[!duplicated(mapping[,1]),]
    names(drugGeneStats) = mapping[, 2]
    gseaRes = fgsea(pathways = pathwaysObj, drugGeneStats, nperm = nperm, nproc = coresNum)
  }else{
    if(grepl("ENSG", pathwaysObj[[1]][1]) == TRUE)
    {
      #users pathways are ensemble ids so drugGeneStats is fine
      gseaRes = fgsea(pathways = pathwaysObj, drugGeneStats, nperm = nperm, nproc = coresNum)
    }else{
      #users pathways are entrez ids so convert drugGeneStats names
      mapping = AnnotationDbi::select(org.Hs.eg.db, keytype="ENSEMBL", keys=as.character(names(drugGeneStats)), columns=c("ENTREZID"))
      mapping = mapping[!duplicated(mapping[,1]),]
      names(drugGeneStats) = mapping[, 2]
      gseaRes = fgsea(pathways = pathwaysObj, drugGeneStats, nperm = nperm, nproc = coresNum)
    }
  }

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
#setwd("C:\\Users\\micha\\Documents\\PMH Research\\GWC CMAP Package")
#load("cmap_sig_rna.RData")
#gseaResDrug = runFgseaPharmDrug("cc", "vinburnine", drug.perturbation, useTstat = TRUE, nperm = 50000)

