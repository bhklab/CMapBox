#' Function to determine which genes drove a drugs connectivity score and create volcano plots that shows which genes drove a drugs connectivity score
#'
#' This function displays volcano plots for the data and drug supplied with the top N and bottom N genes that drove the connectivity score labelled. Top N and bottom N here refer to genes that had positive and negative gene estimates in the data, respectively, and drove the connectivity score the most.
#' Note that which genes drive the connectivity score can vary with the method used to calculate the score. However, in general these genes are the top positive and negative genes that the drug lowered and raised the expression of, respectively.
#' @param genesEst a numeric vector with the estimates (typically LogFC) for the genes
#' @param genesSymb a character vector with the gene symbols for the gene estimates provided
#' @param genesP a numeric vector with the p values for the gene estimates provided
#' @param genesId a character vector with gene Ids for the gene estimates provided. Must be the same type of ID as those used for the row names of the drugPert variable
#' @param drugPert the drug perturbation signature used in the analysis
#' @param drugName the name of the drug in the drug perturbation signature that one would like to inspect
#' @param drugScoreMeth a string specifying which drug repurposing technique was used to score the drugs during each trial. The options for this parameter are currently "gwc" and "fgsea"
#' @param gwcMethod a character string specifying which method was used when computing correlations in the gwc function. The options are spearman or pearson.
#' @param droveNegScore a boolean specifying whether the plots should label the genes that drove the negative connectivity score (reversed the phenotype) between the drug and the data. If FALSE, the genes that drove the positive connectivity score will be labelled.
#' @param drugEst a boolean specifying whether the estimates for each gene of the drug perturbation signature were used in the gwc calculation (TRUE) or if the t-stats for each gene in the drug perturbation signature were used in the gwc calculation (FALSE).
#' @param numbTopGenes number of top positive and negative genes to return information for and include in the plots
#' @param makePlot a boolean specifying whether or not to display the plots showing the genes that the drug interacts with to drive the score. Default is TRUE
#' @param drugVolcPlot a boolean (default is FALSE) specifying whether to generate a volcano plot for the top and bottom 5 genes (by estimate) in the drug that were present in the data supplied. Note that if TRUE the other plots are not generated.
#' @param supressPlot a boolean specifying whether to generate volcano plots that illustrate the genes driving the connectivity score. TRUE by default
#' @return The information used to generate the plots. Specifically, 2 lists.  One that contains the experiment/user supplied data, a subset of this data for the top positive (according to genesEst)
#' genes that drove the score by being altered by the drug, and the same for the most negative genes from the data. The second list is the same but for the drugs data as given by drugPert.
#' The genes will be the same but the info now pertains to how well the drug reduces the expression of these genes.
#' @keywords gwcCMap
#' @export
#' @examples
#' data("geneDataGwc")
#' data("drugPertEx")
#' geneDataClean = mapGenes(geneIds = geneDataGwc$symbol, geneEsts = geneDataGwc$logFC, pvals = geneDataGwc$P.Value, forRankAndPlot = TRUE)
#' drugGeneInfo = plotDrugGeneInteraction(genesEst = geneDataClean$geneEsts, genesSymb = geneDataClean$symbol, genesP = geneDataClean$pvals, genesId = geneDataClean$ensemble, drugPert = drugPertEx, drugName = "BRD-K78431006", drugScoreMeth = "gwcCmapBox", gwcMethod = "pearson", droveNegScore = TRUE, drugEst = TRUE, drugVolcPlot = FALSE, supressPlot = FALSE)
#'

plotDrugGeneInteraction = function(genesEst, genesSymb, genesP, genesId, drugPert, drugName, drugScoreMeth, gwcMethod, droveNegScore, drugEst, numbTopGenes = 5, makePlot = TRUE, drugVolcPlot = FALSE, supressPlot = TRUE)
{
  #  driveGenesNegScore = getSigGeneVolcPlot(genesEst, genesSymb, genesP, genesId, drugPert, topDrugName, gwcMethod, droveNegScore = TRUE, drugEst, drugVolcPlot = FALSE)
  #numbTopGenes = 5
  #genesSymb = geneDataDf$symbol
  #genesP = geneDataDf$pvals
  #genesEst = volcPlotEsts
  #genesId = geneDataDf$ensemble
  #drugName = topDrugName
  #gwcMethod = "pearson"
  #droveNegScore = TRUE
  #drugScoreMeth = "gwcCmapBox"
  #drugName = topDrugName

  volcFrame  = as.data.frame(cbind(genesSymb, genesEst, genesP), stringsAsFactors = FALSE)
  rownames(volcFrame) = genesId
  sameRows = intersect(rownames(drugPert), genesId)
  volcFrame = volcFrame[row.names(volcFrame) %in% sameRows, ]

  colnames(volcFrame) = c("geneSymb", "geneEst", "pvalue")
  volcFrame = transform(volcFrame, geneEst = as.numeric(geneEst))
  volcFrame = transform(volcFrame, pvalue = as.numeric(pvalue))


  drugCol = which(colnames(drugPert) == drugName)
  drugDat = as.data.frame(drugPert[,drugCol ,])
  drugDat = drugDat[row.names(drugDat) %in% sameRows, ]
  drugDat = cbind(drugDat, volcFrame[, "geneSymb"])
  colnames(drugDat)[5] = "geneSymb"


  if(drugScoreMeth == "gwc"){

    truncate.p = 1e-16
    p1 = volcFrame$pvalue
    p1[!is.na(p1) & p1 < truncate.p] <- truncate.p
    p1 = -log10(p1)
    p1 = p1/sum(p1, na.rm = TRUE)

    p2 = drugDat$pvalue
    p2[!is.na(p2) & p2 < truncate.p] <- truncate.p
    p2 = -log10(p2)
    p2 = p2/sum(p2, na.rm = TRUE)

    w = p1 + p2

    if(drugScoreMeth == "gwcCmapBox")
      w = p2

    #for top drug, large neg dataT*pos drugT --> low/drove connectivity score
    #highest rank to highest t-stat lowest rank to lowest t-stat
    #-*- = lowR*lowR, -*+ or +*- = lowR*highR, +*+ = highR*highR
    if(droveNegScore == TRUE){
      direc = FALSE
      direcString = "negative"
      sign = 1
    }else{
      direc = TRUE
      direcString = "positive"
      sign = -1
    }
    if(gwcMethod == "pearson")
    {
      negCol = "blue"
      posCol = "red"
      if(drugEst == FALSE)
      {
        sigGenes = w*volcFrame$geneEst*drugDat$tstat
      }else{
        sigGenes = w*volcFrame$geneEst*drugDat$estimate
      }
      names(sigGenes) = volcFrame$geneSymb
      sigGenesNeg = sigGenes[which(volcFrame$geneEst <= 0)]
      sigGenesPos = sigGenes[which(volcFrame$geneEst > 0)]
    }else{
      negCol = "red"
      posCol = "blue"
      if(drugEst == FALSE){
        sigGenesNeg = w*rank(sign*-volcFrame$geneEst)*rank(drugDat$tstat)
        sigGenesPos = w*rank(sign*volcFrame$geneEst)*rank(-drugDat$tstat)
      }else{
        sigGenesNeg = w*rank(sign*-volcFrame$geneEst)*rank(drugDat$estimate)
        sigGenesPos = w*rank(sign*volcFrame$geneEst)*rank(-drugDat$estimate)
      }
      names(sigGenesNeg) = volcFrame$geneSymb
      names(sigGenesPos) = volcFrame$geneSymb
    }
    #choosing lowest sigGenesNeg and sigGenesPos to represent genes that drove -ve connectivity score
    #so, for low rank,

    #dataT < 0 --> these are negative genes of data that with drug gave high score
    #high expression in disease group and drug lowers expression of gene

    #Spear sigGenesNeg = w*rank(-volcFrame$tstat)*rank(topDrugDat$tstat)
    #Spear names(sigGenesNeg) = volcFrame$geneSymb
    #dataT > 0 --> these are +ve genes of data that with drug gave high score
    #low expression in disease group comp normal group and drug raises expression

    #SPEAR sigGenesPos = w*rank(volcFrame$tstat)*rank(-topDrugDat$tstat)
    # SPEAR names(sigGenesPos) = volcFrame$geneSymb

    topNeg = names(sort(sigGenesNeg, index.return = TRUE, decreasing = direc)$x)[1:numbTopGenes]
    topNegInds = c()
    for(i in 1:length(topNeg))
      topNegInds = c(topNegInds, which(volcFrame$geneSymb == topNeg[i]))
    negSubsetDat = volcFrame[topNegInds, ]
    negSubsetDrug = drugDat[topNegInds, ]

    topPos = names(sort(sigGenesPos, index.return = TRUE, decreasing = direc)$x)[1:numbTopGenes]
    topPosInds = c()
    for(i in 1:length(topPos))
      topPosInds = c(topPosInds, which(volcFrame$geneSymb == topPos[i]))
    posSubsetDat = volcFrame[topPosInds, ]
    posSubsetDrug = drugDat[topPosInds, ]

  }else if(drugScoreMeth == "xsum" | drugScoreMeth == "fgsea" | drugScoreMeth == "gwcCmapBox"){
    #does this work for fgsea and gwcCmapBox? have sign only of geneEst and this just takes most + and - drugDat
    #genes with correct signs. likely works
    #fgsea each gene equal impact so free to choose any genes so long as sign is correct
    #gwcCmapBox considers magnitude so just take most neg or  pos, since free to choose for fgsea
    #most sig genes in xsum likely same as most sig genes in gwcCmapBox so works as well
    #put different string in plot, or outside function, depending on method

    if(droveNegScore == TRUE){
      direc = FALSE
      direcString = "negative"
      sign = 1
    }else{
      direc = TRUE
      direcString = "positive"
      sign = -1
    }

    x1 = volcFrame$geneEst
    names(x1) = volcFrame$geneSymb
    if(drugEst == TRUE)
      y1 = drugDat$estimate
    if(drugEst == FALSE)
      y1 = drugDat$tstat
    names(y1) = drugDat$geneSymb
    sumUpX = 0
    sumDownX = 0
    upInDis = x1[which(x1 >= 0)]
    downInDis = x1[which(x1 < 0)]
    topAndBotCompound = c(y1[order(y1)][1:length(upInDis)], y1[order(y1, decreasing = TRUE)][1:length(downInDis)])
    upInDisInter = intersect(names(upInDis), names(topAndBotCompound))
    downInDisInter = intersect(names(downInDis), names(topAndBotCompound))

    if(length(upInDisInter) > 0)
      upX = topAndBotCompound[upInDisInter]
    if(length(downInDisInter) > 0)
      downX = topAndBotCompound[downInDisInter]

    topNeg = names(sort(sign*upX)[1:numbTopGenes])
    topNegInds = c()
    for(i in 1:length(topNeg))
      topNegInds = c(topNegInds, which(volcFrame$geneSymb == topNeg[i]))
    negSubsetDat = volcFrame[topNegInds, ]
    negSubsetDrug = drugDat[topNegInds, ]

    topPos = names(sort(sign*downX, decreasing = TRUE)[1:numbTopGenes])
    topPosInds = c()
    for(i in 1:length(topPos))
      topPosInds = c(topPosInds, which(volcFrame$geneSymb == topPos[i]))
    posSubsetDat = volcFrame[topPosInds, ]
    posSubsetDrug = drugDat[topPosInds, ]

    negCol = "red"
    posCol = "blue"
  }

  if(drugVolcPlot == TRUE)
  {
    mostPosInds = sort(drugDat[,"estimate"], index.return = TRUE, decreasing = TRUE)$ix
    posSubset = drugDat[mostPosInds[1:numbTopGenes], ]
    mostNegInds = sort(drugDat[,"estimate"], index.return = TRUE)$ix
    negSubset = drugDat[mostNegInds[1:numbTopGenes], ]

    with(drugDat, plot(estimate, -log10(pvalue), pch=20, main=paste("Drug (",drugName, ") Volcano Plot"), xlim=c(1.4*min(drugDat$estimate), 1.4*max(drugDat$estimate)), font.axis = 2, cex.lab = 1.3, cex.axis = 1.15))
    with(negSubset, points(estimate, -log10(pvalue), pch=20, col="blue", cex = 2))
    with(negSubset, textxy(estimate, -log10(pvalue), labs = geneSymb, cex = 1.5, font = 2, col = "green4"))
    with(posSubset, points(estimate, -log10(pvalue), pch=20, col="red", cex = 2))
    with(posSubset, textxy(estimate, -log10(pvalue), labs = geneSymb, cex = 1.5, font = 2, col = "green4"))

  }else{

    #pdf(paste("sig gene queried", genesSigStart,"to",genesSigEnd,"by",genesSigInc,"Volcano Plot of Data Top Drug Genes.pdf"))
    #was previously at < 2 probably just an error though
    par(mfrow=c(2,1))
    if(length(unique(genesEst)) > 2){
      with(volcFrame, plot(geneEst, -log10(pvalue), pch=20, main=paste("Volcano Plot for the Data with the genes that drove the", direcString, "connectivity score between drug", drugName, "and the data labelled" ), xlim=c(1.2*min(volcFrame$geneEst), 1.2*max(volcFrame$geneEst)), font.axis = 2, cex.lab = 1.3, cex.axis = 1.15))
      with(negSubsetDat, points(geneEst, -log10(pvalue), pch=20, col=negCol, cex = 2))
      with(negSubsetDat, textxy(geneEst, -log10(pvalue), labs = geneSymb, cex = 1.5, font = 2, col = "green4"))
      with(posSubsetDat, points(geneEst, -log10(pvalue), pch=20, col=posCol, cex = 2))
      with(posSubsetDat, textxy(geneEst, -log10(pvalue), labs = geneSymb, cex = 1.5, font = 2, col = "green4"))
      legend(bty = "n", 0.8*max(volcFrame$geneEst), 0.97*max(-log10(volcFrame$pvalue)), legend=c("Low in Disease", "High in Disease"), col=c("blue", "red"), pch = 19, cex=1.1, pt.cex = 1.3, x.intersp = .75)
      #dev.off()
    }else{
      mostInfluentialGenes = c(1:numbTopGenes)
      colnames(negSubsetDat)[2] = c("gene direction in signature")
      colnames(posSubsetDat)[2] = c("gene direction in signature")
      with(negSubsetDat, plot(`gene direction in signature`, mostInfluentialGenes, pch=20, main=paste("Plot for the Data with the genes that drove the", direcString, "connectivity score between drug", drugName, "and the data labelled" ), xlim=c(1.2*min(volcFrame$geneEst), 1.2*max(volcFrame$geneEst)), font.axis = 2, cex.lab = 1.3, cex.axis = 1.15))
      with(negSubsetDat, points(`gene direction in signature`, mostInfluentialGenes, pch=20, col=negCol, cex = 2))
      with(negSubsetDat, textxy(`gene direction in signature`, mostInfluentialGenes, labs = geneSymb, cex = 1.5, font = 2, col = "green4"))
      with(posSubsetDat, points(`gene direction in signature`, mostInfluentialGenes, pch=20, col=posCol, cex = 2))
      with(posSubsetDat, textxy(`gene direction in signature`, mostInfluentialGenes, labs = geneSymb, cex = 1.5, font = 2, col = "green4"))
      legend(bty = "n", 0, 3, legend=c("Low in Disease", "High in Disease"), col=c("blue", "red"), pch = 19, cex=1.1, pt.cex = 1.3, x.intersp = .75)
    }

    #pdf(paste("sig gene queried", genesSigStart,"to",genesSigEnd,"by",genesSigInc,"Volcano Plot Top Drug Narrow X.pdf"))
    with(drugDat, plot(estimate, -log10(pvalue), pch=20, main=paste("Drug (",drugName, ") Volcano Plot with the genes that drove the", direcString, "connectivity score bewteen the drug and the data labelled"), xlim=c(1.4*min(drugDat$estimate), 1.4*max(drugDat$estimate)), font.axis = 2, cex.lab = 1.3, cex.axis = 1.15))
    with(negSubsetDrug, points(estimate, -log10(pvalue), pch=20, col=negCol, cex = 2))
    with(negSubsetDrug, textxy(estimate, -log10(pvalue), labs = geneSymb, cex = 1.5, font = 2, col = "green4"))
    with(posSubsetDrug, points(estimate, -log10(pvalue), pch=20, col=posCol, cex = 2))
    with(posSubsetDrug, textxy(estimate, -log10(pvalue), labs = geneSymb, cex = 1.5, font = 2, col = "green4"))
    legend(bty = "n", 0.8*max(drugDat$estimate), 0.97*max(-log10(drugDat$pvalue)), legend=c("Low in Disease", "High in Disease"), col=c("blue", "red"), pch = 19, cex=1.1, pt.cex = 1.3, x.intersp = .75)
    #dev.off()
    par(mfrow=c(1,1))


  }
  genesDroveScore = c(topPos, topNeg)
  experimentList = list(volcFrame, posSubsetDat, negSubsetDat)
  drugList = list(drugDat, posSubsetDrug, negSubsetDrug)
  names(experimentList) = c("experimentData", "posGeneDriversExper", "negGeneDriversExper")
  names(drugList) = c("drugData", "posGeneDriversDrug", "negGeneDriversDrug")
  retList = list(experimentList, drugList)
  names(retList)= c("experimentList", "drugList")

  return(retList)

}
