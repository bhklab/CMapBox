#' Function to generate heatmaps that show a drugs influence on a set of genes
#'
#' This function returns a heat map that shows how a drug changes the expression of a set of genes, as determined by molecular profiles in the pharmacoSet supplied.
#' @param drugPert the drug perturbation signature of interest
#' @param pharmSet the pharmacoSet used to generate the drug perturbation signature
#' @param mDataType the type of molecular data to obtain from the PharmacoSet (see PharmacoGx::molecularProfiles)
#' @param drugName the name of the drug one would like to inspect for its effect on the supplied genes
#' @param genesOfInt a character vector containing the names of the genes one would like to observe the drugs effect on
#' @param makePlot a boolean specifying whether to generate the plot showing the gene expression changes due to the drug. Default is TRUE
#' @return a list containing 1. beforeDrugMat which is a matrix containing the expression levels for the genes before the drug was applied.
#'  2. afterDrugMat which is a matrix containing the expression levels for the genes before the drug was applied. 3. geneChangePlot, which is a
#'  levelplot object showing the change in expression for the genes due to the drugs
#' @export
#' @examples
#' data("drugPertEx")
#' data("psetSub")
#' #example in getSigGeneVolPlot function yields below genes
#' drivingGenes = c("ABCC5", "VPS28", "RPS6KA2", "TRIM2", "CD47", "MELK", "CDKN3", "GGH", "MCM7", "KIF2C")
#' plotAndMats = getDrugsImpactOnGenes(drugPert = drugPertEx, pharmSet = psetSub, mDataType = "rna", drugName = "BRD-K78431006", genesOfInt = drivingGenes)
#'

getDrugsImpactOnGenes = function(drugPert, pharmSet, mDataType, drugName, genesOfInt, makePlot = TRUE)
{
  #pharmSet = L1000_compounds
  #mDataType = "rna"
  #drugName = "BRD-K19687926"
  #genesOfInt = c("MCOLN1" ,"JUN"  ,  "CERK"  , "IGFBP3", "TIMP2" , "E2F2" ,  "KIF14" , "PLK1",   "DCK" ,   "UBE2C" )

  #drugName = topDrugName
  #genesOfInt = driveGenesNegScore

  drugCol = which(colnames(drugPert) == drugName)
  drugDat = as.data.frame(drugPert[,drugCol ,])
  geneEns = rownames(drugPert)
  mapping = AnnotationDbi::select(org.Hs.eg.db, keytype="ENSEMBL", keys=as.character(geneEns), columns=c("SYMBOL"))
  mapping = mapping[!duplicated(mapping[,1]),]
  drugDat[, ncol(drugDat) + 1] = mapping[, 2]
  colnames(drugDat)[ncol(drugDat)] = "geneSymb"

  inds = c()
  for(i in 1:length(genesOfInt))
    inds = c(inds, which(drugDat$geneSymb == genesOfInt[i]))
  subsetDrug = drugDat[inds, ]

  pharmDat = molecularProfiles(pharmSet, mDataType)
  #note cant change rownames(pharmDat) due to memory issues
  if(sum(grepl("_at", rownames(pharmDat))) > nrow(pharmDat)/2)
    rownamesPharmDat = gsub("_at_", "", paste0(rownames(pharmDat), "_"))
  if(sum(grepl("geneid.", rownames(pharmDat))) > nrow(pharmDat)/2)
  {
    geneEntrezs = gsub("geneid.", "", rownames(pharmDat))
    mapping = AnnotationDbi::select(org.Hs.eg.db, keytype="ENTREZID", keys=as.character(geneEntrezs), columns=c("ENSEMBL"))
    mapping = mapping[!duplicated(mapping[,1]),]
    rownamesPharmDat = mapping[, 2]
  }
  colsInfo = phenoInfo(pharmSet, mDataType)
  contRows = which(colsInfo$xptype == "control")
  contDat = colsInfo[contRows, ]

  drugColsInfo = colsInfo[which(colsInfo$drugid == drugName), ]
  drugBatchIds = drugColsInfo$batchid
  drugContRows = which(contDat$batchid %in% drugBatchIds)

  contMat = as.data.frame(NULL)
  for(i in 1:length(drugBatchIds))
  {
    contBatchId = drugBatchIds[i]
    contRows = which(contDat$batchid %in% contBatchId)
    nameVec = as.character(rownames(contDat)[contRows])
    #nameVec = as.character(contDat$samplename[contRows])
    exprDat = 0
    for(j in 1:length(nameVec))
    {
      exprDat = exprDat + pharmDat[, which(colnames(pharmDat) == nameVec[j])]/length(nameVec)
    }
    contMat = rbind(contMat, exprDat)
  }
  #rownames(contMat) = drugBatchIds
  if(ncol(contMat) > 1){
    colnames(contMat) = rownamesPharmDat
  }else{
    names(contMat) = rownamesPharmDat
  }
  contMat = t(contMat)
  drugMat = pharmDat[, which(colnames(pharmDat) %in% rownames(drugColsInfo))]
  if(!is.null(dim(drugMat))){
    colnames(drugMat) = drugBatchIds
  }else{
    names(drugMat) = drugBatchIds
  }

  ensIds = c(rownames(subsetDrug))
  geneNames = as.character(subsetDrug$geneSymb)
  geneRows = which(rownames(contMat) %in% ensIds)
  #just use first batch results
  contMatSub = contMat[geneRows, ]

  if(!is.null(dim(drugMat))){
    drugMatSub = drugMat[geneRows, ]
    rownames(contMatSub) = geneNames
    rownames(drugMatSub) = geneNames
    #need to have labels and title in d3heatmap before using it
    #plot(d3heatmap(t((drugMatSub/contMatSub - 1)*100), scale="column", colors=colorRampPalette(c("blue", "red"))( 4 ), dendrogram = "none"), xlab = "Batch ID")
    #dev.new()
    par(mfrow=c(1,2))
    plot(levelplot(t(contMatSub), scales=list(x=list(rot=90)), cex.lab = 1.6, xlab = "Control Batch ID", ylab = "Genes", main = paste("Gene Expression levels before drug", drugName, "was Applied According to the pharmacoSet data"), aspect = "fill"))
    beforeDrugPlot = levelplot(t(contMatSub), scales=list(x=list(rot=90)), cex.lab = 1.6, xlab = "Control Batch ID", ylab = "Genes", main = paste("Gene Expression levels before drug", drugName, "was Applied According to the pharmacoSet data"), aspect = "fill")
    beforeDrugMat = t(contMatSub)
    #dev.new()
    plot(levelplot(t(drugMatSub), scales=list(x=list(rot=90)), cex.lab = 1.6, xlab = "Batch ID", ylab = "Genes", main = paste("Gene Expression levels after drug", drugName, "was Applied According to the pharmacoSet data"), aspect = "fill"))
    afterDrugPlot = levelplot(t(drugMatSub), scales=list(x=list(rot=90)), cex.lab = 1.6, xlab = "Batch ID", ylab = "Genes", main = paste("Gene Expression levels after drug", drugName, "was Applied According to the pharmacoSet data"), aspect = "fill")
    afterDrugMat = t(drugMatSub)
    frameToPlot = t((drugMatSub/contMatSub - 1)*100)
    legRange = max(abs(frameToPlot))
    #dev.new()
    if(makePlot == TRUE) plot(levelplot(frameToPlot, at=seq(-legRange, legRange, length.out = 100), scales=list(x=list(rot=90)), cex.lab = 1.6, xlab = "Batch ID", ylab = "Genes", main = paste("Percent change in Gene Expression due to drug", drugName, "in multiple batches According to the pharmacoSet data"), aspect = "fill"))
    geneChangePlot = levelplot(frameToPlot, at=seq(-legRange, legRange, length.out = 100), scales=list(x=list(rot=90)), cex.lab = 1.6, xlab = "Batch ID", ylab = "Genes", main = paste("Percent change in Gene Expression due to drug", drugName, "in multiple batches According to the pharmacoSet data"), aspect = "fill")
  }else{
    drugMatSub = drugMat[geneRows]
    names(contMatSub) = geneNames
    names(drugMatSub) = geneNames
    before_drug = contMatSub
    after_drug = drugMatSub
    #plot(levelplot(t(before_drug), scales=list(x=list(rot=90)), cex.lab = 1.6, xlab = "Control Batch ID", ylab = "Genes", main = paste("Gene Expression levels before drug", drugName, "was Applied According to the pharmacoSet data"), aspect = "fill"))
    #plot(levelplot(t(after_drug), scales=list(x=list(rot=90)), cex.lab = 1.6, xlab = "Batch ID", ylab = "Genes", main = paste("Gene Expression levels after drug", drugName, "was Applied According to the pharmacoSet data"), aspect = "fill"))
    frameToPlot = t(cbind(before_drug, after_drug))
    legRange = max(abs(frameToPlot))
    #dev.new()
    if(makePlot == TRUE) plot(levelplot(frameToPlot, at=seq(-legRange, legRange, length.out = 100), scales=list(x=list(rot=90)), cex.lab = 1.6, xlab = "Expression Levels Before the Drug (Far Left) and after the Drug was Applied", ylab = "Genes", main = paste("Change in Gene Expression due to drug", drugName, "According to the pharmacoSet data (note: only 1 batch found)"), aspect = "fill"))
    geneChangePlot = levelplot(frameToPlot, at=seq(-legRange, legRange, length.out = 100), scales=list(x=list(rot=90)), cex.lab = 1.6, xlab = "Expression Levels Before the Drug (Far Left) and after the Drug was Applied", ylab = "Genes", main = paste("Change in Gene Expression due to drug", drugName, "According to the pharmacoSet data (note: only 1 batch found)"), aspect = "fill")
    beforeDrugMat = contMatSub
    afterDrugMat = drugMatSub

  }

  retList = list(beforeDrugMat, afterDrugMat, geneChangePlot)
  names(retList) = c("beforeDrugMat", "afterDrugMat", "geneChangePlot")
  return(retList)
}

#cellExpFrame = t(as.data.frame(pertNumber(pharmSet)))
#rownames(cellExpFrame) = drugInfo(pharmSet)$pert_iname
