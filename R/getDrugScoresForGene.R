
#' Function to use specified genes in a data set to get connectivity scores and p values for drugs via a genome wide correlation
#'
#' This function will use the genes inputed to determine connectivity scores and p values for the drugs in a drug perturbation signature
#' @param drugNameVec a character vector containing the names of the drugs in the drug perturbation signature to compute a more accurate p value for
#' @param geneIds a character vector containing the gene symbols, ensemble IDs, or entrez IDs for the genes to be analyzed
#' @param geneEsts a numeric vector containing estimates for the difference in expression between the two phenotyes. Note that positive values of this estimate (tstat, logFC, etc) must correspond to higher gene expression in the phenotype one would like to reverse
#' @param pvals P values from a t-test assessing the diferential expression of the genes between the two phenotypes
#' @param drugPert a drug perturbation signature (object of class PharmacoSig) with rownames that correspond to the ensemble IDs of the genes
#' @param pharmSet the PharmacoSet used to generate the drug perturbation signature. Supplying this will add additonal info about the drugs in the results table.
#' @param drugScoreMeth a string specifying which drug repurposing technique to use to score the drugs during each trial. The options for this parameter are currently "gwc" and "fgsea". Default is "gwc"
#' @param numbPerms The number of permutations to be used to compute the p value in the gwc function. 
#' @param gwcMethod a character string specifying which method to use when computing correlations in the gwc function. The options are spearman (default) or pearson.
#' @param drugEst a boolean specifying whether to use the estimates for each gene of the drug perturbation signature in the gwc calculation (TRUE) or to use the t-stats for each gene in the drug perturbation signatur ein the gwc calculation (FALSE). Default is TRUE
#' @return a data frame containing scores and p values for the drugs when analyzed using the supplied genes
#' @keywords gwcCMap
#' @export
#' @examples
#' data("psetSub")
#' data("geneDataGwc")
#' data("drugPertEx")
#' geneDataClean = cleanData(geneIds = geneDataGwc$ensembl_id, geneEsts = geneDataGwc$logFC, pvals = geneDataGwc$P.Value)
#' #keep genes with p value < 0.10 and |logFC| > 0.5
#' rowsToKeep = removeGenes(idsOne = geneDataClean$ensemble, idsTwo = geneDataClean$ensemble, valuesTwo = cbind(geneDataClean$pvals, abs(geneDataClean$geneEsts)), conditionVals = c(0.10, .5), conditionsDirec = c(FALSE, TRUE))
#' sigGeneData = geneDataClean[rowsToKeep, ]
#' gwcResults = getDrugScoresForGenes(geneIds = sigGeneData$ensemble, geneEsts = sigGeneData$geneEsts, pvals = sigGeneData$pvals, drugPert = drugPertEx, pharmSet = psetSub)
#' 

getDrugScoresForGenes = function(geneIds, geneEsts, pvals, drugPert = NULL, pharmSet = NULL, drugScoreMeth = "gwc", gwcMethod = "spearman", numbPerms = 10000, drugEst = TRUE)
{
  #convert geneEsts and pvals in case they are factors or character vectors
  geneEsts = as.numeric(as.character(geneEsts))
  pvals = as.numeric(as.character(pvals))
  
  geneDataDf = cleanData(geneIds, geneEsts, pvals, forRankAndPlot = TRUE)
  
  #pharmacoGx doesnt appear to install if needed, below installs it if it hasn't been before
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
  
  genesCmap = genesCmapAll
  ordering = order(abs(genesCmap$pvalue), decreasing = FALSE)
  
  drugNotes = NULL
  if(!is.null(pharmSet))
    drugNotes = drugInfo(pharmSet)
  
    genesCmapUse = genesCmapAll
    
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
    for(i in 1:numDrugs)
    {
      if(drugEst == TRUE)
      {
        drugsList[[i]] = drugPert[,i ,c("tstat", "pvalue")]
      }else{
        drugsList[[i]] = drugPert[,i ,c("estimate", "pvalue")]
      }
    }
    
    csFunc = function(z) connectivityScore(x = genesCmapUse, y = z,method=drugScoreMeth, gwc.method=gwcMethod, nperm = numbPerms)
    environment(csFunc) <- .GlobalEnv
    
    clusterExport(cl, list("drugsList", "genesCmapUse", "numbPerms", "gwcMethod", "drugScoreMeth", "connectivityScore", "gwc", "intersectList", "pco", "corWeighted", "combineTest"), envir = environment())
    dataCor = clusterApplyLB(cl, drugsList, csFunc)
    stopCluster(cl)
    
    for(i in 1:numDrugs)
    {
      if(!is.null(drugNotes))
      {
        drugInd = which(rownames(drugNotes) == drugNames[i])
        drugInfoTab = drugNotes[drugInd, ]
        cmapTab = rbind(cmapTab, c(dataCor[[i]], as.character(drugInfoTab)), stringsAsFactors = FALSE) 
        colnames(cmapTab) = c("gwc score", "fdr adjusted pvalue", colnames(drugNotes[1,]))
      }else{
        cmapTab = rbind(cmapTab, c(dataCor[[i]]), stringsAsFactors = FALSE) 
        colnames(cmapTab) = c("gwc score", "fdr adjusted pvalue")
      }
    }
    
    rownames(cmapTab) = drugNames[1:numDrugs]
    sorting = order(cmapTab[,"gwc score"], decreasing = FALSE)
    cmapTabSort = cmapTab[sorting,]
    cmapTabSort[,2] = p.adjust(cmapTabSort[,2], method = "fdr")
    
    #one of the drugs score kept showing up as NA, so remove it
    naDrugRows = which(is.na(cmapTabSort[,1]))
    if(length(naDrugRows) > 0)
      cmapTabSort = cmapTabSort[-naDrugRows, ]
    
  return(cmapTabSort)
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


