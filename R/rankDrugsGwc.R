
#' Function to identify drug candidates.
#'
#' This function ranks a list of drugs according to their ability to reverse a phenotype. This is done by first computing scores for the drugs according to one of the various drug repurposing techniques over multiple trials, where each trial
#' contains a different number of differentially expressed genes between the disease and normal phenotype. At each trial, the drugs are ranked relative to one another in terms of their ability to reverse the disease phenotype. Additionally,
#' drugs are ranked according to how stable they are from trial to trial (as one varies the number of genes). A final score is given to each drug and is determined according to a drugs rank throughout each trial and its change in rank from trial to trial. For N drugs,
#' a final score of -N would imply the drug was the best at reversing the disease phenotype every single trial and that the drugs rank changed the least from one trial to the next, thus meaning it is both the best drug at reversing the disease phenotype and the most stable drug candidate.
#' @param geneIds a character vector containing the gene symbols, ensemble IDs, or entrez IDs for the genes to be analyzed
#' @param geneEsts a numeric vector containing estimates or directions (+1 -> up or -1 -> down) for the difference in expression between the two phenotyes. Note that positive values of this estimate (+1, or tstat, logFC, etc) must correspond to higher gene expression in the phenotype one would like to reverse
#' @param pvals P values from a t-test assessing the diferential expression of the genes between the two phenotypes. If not provided (i.e direction information only for geneEsts) then volcano plots for the data will not be generated and repeat gene IDs will be removed arbitrarily (first ID kept) as opposed to by p value
#' @param drugPert a drug perturbation signature (object of class PharmacoSig) with rownames that correspond to the ensemble IDs of the genes
#' @param pharmSet the PharmacoSet used to generate the drug perturbation signature. Supplying this will add additonal info about the drugs in the ranking tables.
#' @param drugScoreMeth a string specifying which drug repurposing technique to use to score the drugs during each trial. The options for this parameter are currently "gwc", "gwcCmapBox", "fgsea", and "xsum". Default is "gwcCmapBox"
#' @param pCut a boolean specifying whether to use pvalues to remove insignificant genes from the CMAP analysis (TRUE) or to remove genes from the CMAP analysis according to their supplied gene estimates (FALSE). Initially TRUE, but Defaults to FALSE if no p-values are supplied
#' @param cutOff if pCut is TRUE then this value represents the p value threshold used to filter out genes. If pCut is FALSE then this value represents the fraction of genes present in both the data and drug perturbation signature with the top absolute value of gene estimates that will be left in the analysis. cutOff should be between 0 and 1. If p values ar enot provided, all genes supplied will be used in the analysis
#' @param genesToStart a value between 0 and 1 representing the fraction of genes present in the data and drug perturbation signature to use in the first iteration of the CMAP analysis. Recommended to be at least 0.40 to avoid large changes during the early iterations due to having too few genes present.
#' @param numbIters The number of iterations that will occur in the analysis. numbIters will set the rate at which genes are added to the analysis.
#' @param numbPerms The number of permutations to be used to compute the p value in the drug repurposing functions.
#' @param gwcMethod a character string specifying which method to use when computing correlations in the gwc function. The options are spearman (default) or pearson.
#' @param volcPlotEsts a numeric vector containing gene estimates to be used for volcano plots, if not provided geneEsts will be used. If one desires to use t-stats for geneEsts in the gwc analysis it is recommended that one supply the logFC of the genes here if volcano plots are desired.
#' @param drugEst a boolean specifying whether to use the estimates for each gene of the drug perturbation signature in the calculations (TRUE) or to use the t-stats for each gene in the drug perturbation signatur ein the calculations (FALSE). Default is TRUE
#' @param inspectResults a boolean specifying whether to display plots that will allow one to check that the correct drugs have been selected based on the data supplied. Default is TRUE (show plots)
#' @param showMimic a boolean (default is FALSE) specifying whether to show plots for the drug that upregulates and down regulates the genes that are overexpressed and under expressed in the disease state.
#' @param mDataType a string specifying the type of molecular data to retrieve molecular profiles for from the pharmacoSet, if one desires to inspect the results of the analysis (inspectResults = TRUE). Default is "rna"
#' @param extraData a data.frame where each column represents values one would like to inspect to determine if the genes corresponding to these values should be removed if they meet the conditions specified in extraCut and extraDirec. Useful if one would like to remove genes based on logFC or other values. extra data must have the same number of rows as the length of vectors geneIds, geneEsts, and pvals.
#' @param extraCut a numeric vector, with the nth value corresponding to the nth column of the m by n data.frame supplied in extraData, where each value indicates the value that needs to be reached for the columns of the data fram in order for the ids with that value to be kept or removed. Whether removal occurs when a value is greater or less than the value in this vector is specified in the extraDirec variable.
#' @param extraDirec a boolean vector, where TRUE means that ids whose value is greater than that specified in extraVals will be removed and FALSE means those with a value less than conditionVals will be removed.
#' @param drugNameVec a character vector specifying the names of the drugs to run through the pipeline. Useful for saving time by running a higher nperm analysis to get more reliable p values on the top drugs identified in a smaller nperm analysis with all the drugs. Default is NULL and all drugs are tested.
#' @return a data frame with information about the drugs in the analysis, including the drugs final scores. There are multiple connectivity scores, p values, and fdr adjusted p values in the results table as each value corresponds to one iteration in the analysis.
#' @keywords gwcCMap
#' @import Biobase org.Hs.eg.db piano parallel lattice calibrate AnnotationDbi stats PharmacoGx
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#' @importFrom piano runGSA
#' @importFrom piano loadGSC
#' @importFrom stats complete.cases
#' @importFrom parallel detectCores makeCluster clusterExport clusterApplyLB stopCluster
#' @importFrom PharmacoGx gwc molecularProfiles phenoInfo intersectList drugInfo
#' @importFrom lattice levelplot
#' @importFrom calibrate textxy
#' @importFrom Biobase note
#' @export
#' @examples
#' data("geneDataGwc")
#' data("drugPertEx")
#'
#' data("psetSub")
#' #use below command to use all of the CMAP drugs in your analysis
#' #library(PharmacoGx)
#' #drugPertEx = downloadPSet("CMAP")
#' drugResults = rankDrugsGwc(geneIds = geneDataGwc$symbol, geneEsts = geneDataGwc$t, pvals = geneDataGwc$P.Value, drugPert = drugPertEx, pharmSet = psetSub, drugScoreMeth = "gwc", volcPlotEsts = geneDataGwc$logFC)
#' #below line shows how to additionally filter genes based on logFC
#' #drugResults = rankDrugsGwc(geneIds = geneDataGwc$ensembl_id, geneEsts = geneDataGwc$t, pvals = geneDataGwc$P.Value, drugPert = drugPertEx, pharmSet = psetSub, volcPlotEsts = geneDataGwc$logFC, extraData = cbind(abs(geneDataGwc$logFC)), extraCut = c(0.5), extraDirec = c(TRUE))

rankDrugsGwc = function(geneIds, geneEsts, drugPert, pharmSet = NULL, pvals = NULL, drugScoreMeth = "gwcCmapBox", pCut = TRUE, cutOff = 0.05, genesToStart = 0.20, numbIters = 10, gwcMethod = "spearman", numbPerms = 1000, volcPlotEsts = NA, drugEst = TRUE, inspectResults = FALSE, showMimic = FALSE, mDataType = "rna", extraData = NA, extraCut = NA, extraDirec = NA, drugNameVec = NULL)
{
  #########################################Section 1: Data Preperation#######################################################
  #To Do
  #1. add pharmacoGx to imports and remove combineTest from here and remove corWeighted from gwcCmapBox + here and see if this still runs
  #2. look at old package getDrugPVals and add method/option of getting p values in this func for top drugs only so its quick.
  #   may be as simple as adding option to give drug names into function so only use those drugs on second run through with higher nperm
  #3. make sure drug ran info is returned somehow, may already be in data framr eturned if recalling corrrectly

  if(length(geneIds) != length(unique(geneIds)))
     stop("There are duplicated or missing geneIds. Please ensure there are no NA ids and that there is only 1 of each gene in the provided data as we do not know which of the repeat genes should be
          used in the analysis")
  if(1)
  {

  combineTest <-
    function(p, weight, method=c("fisher", "z.transform", "logit"), hetero=FALSE, na.rm=FALSE) {
      if(hetero) { stop("function to deal with heterogeneity is not implemented yet!") }
      method <- match.arg(method)
      na.ix <- is.na(p)
      if(any(na.ix) && !na.rm) { stop("missing values are present!") }
      if(all(na.ix)) { return(NA) } ## all p-values are missing
      p <- p[!na.ix]
      k <- length(p)
      if(k == 1) { return(p) }
      if(missing(weight)) { weight <- rep(1, k); }
      switch(method,
             "fisher"={
               cp <- pchisq(-2 * sum(log(p)), df=2*k, lower.tail=FALSE)
             },
             "z.transform"={
               z <- qnorm(p, lower.tail=FALSE)
               cp <- pnorm(sum(weight * z) / sqrt(sum(weight^2)), lower.tail=FALSE)
             },
             "logit"={
               tt <- (- sum(log(p / (1 - p)))) / sqrt(k * pi^2 * (5 * k + 2) / (3 * (5 * k + 4)))
               cp <- pt(tt,df=5*k+4, lower.tail=FALSE)
             })
      return(cp)
    }

  corWeighted <-
    function (x, y, w, method=c("pearson", "spearman"), alternative=c("two.sided", "greater", "less"), nperm=0, nthread=1, setseed, na.rm=FALSE) {

      ######################
      wcor <- function (d, w, na.rm=TRUE) {
        ### NOTE::: THIS FORMULA CAN SUFFER CATASTROPHIC CANCELATION AND SHOULD BE FIXED!!!
        #     s <- sum(w, na.rm=na.rm)
        #     m1 <- sum(d[ , 1L] * w, na.rm=na.rm) / s
        #     m2 <- sum(d[ , 2L] * w, na.rm=na.rm) / s
        #     res <- (sum(d[ , 1L] * d[ , 2L] * w, na.rm=na.rm) / s - m1 * m2) / sqrt((sum(d[ , 1L]^2 * w, na.rm=na.rm) / s - m1^2) * (sum(d[ , 2L]^2 * w, na.rm=na.rm) / s - m2^2))
        CovM <- cov.wt(d, wt=w)[["cov"]]
        res <- CovM[1,2]/sqrt(CovM[1,1]*CovM[2,2])
        return (res)
      }

      ######################

      if (missing(w)) { w <- rep(1, length(x)) / length(x) }
      if (length(x) != length(y) || length(x) != length(w)) { stop("x, y, and w must have the same length") }
      method <- match.arg(method)
      if (method == "spearman") {
        x <- rank(x)
        y <- rank(y)
      }
      alternative <- match.arg(alternative)

      res <- c("rho"=NA, "p"=NA)

      ## remove missing values
      ccix <- complete.cases(x, y, w)
      if(!all(ccix) && !na.rm) { warning("Missing values are present") }
      if(sum(ccix) < 3) {
        return(res)
      }
      x <- x[ccix]
      y <- y[ccix]
      w <- w[ccix]

      wc <- wcor(d=cbind(x, y), w=w)
      res["rho"] <- wc
      if (nperm > 1) {
        #if (!missing(setseed)) { set.seed(setseed) }
        splitix <- parallel::splitIndices(nx=nperm, ncl=nthread)
        if (!is.list(splitix)) { splitix <- list(splitix) }
        splitix <- splitix[sapply(splitix, length) > 0]
        mcres <- parallel::mclapply(splitix, function(x, xx, yy, ww) {
          pres <- sapply(x, function(x, xx, yy, ww) {
            ## permute the data and the weights
            d2 <- cbind(xx[sample(1:length(xx))], yy[sample(1:length(yy))])
            w2 <- ww[sample(1:length(ww))]
            return(wcor(d=d2, w=w2))
          }, xx=xx, yy=yy, ww=ww)
          return(pres)
        }, xx=x, yy=y, ww=w)
        perms <- do.call(c, mcres)

        switch (alternative,
                "two.sided" = {
                  if (res["rho"] < 0) { p <- sum(perms <= res, na.rm=TRUE) } else { p <- sum(perms >= res, na.rm=TRUE) }
                  if (p == 0) { p <- 1 / (nperm + 1) } else { p <- p / nperm }
                  p <- p * 2
                },
                "greater" = {
                  p <- sum(perms >= res, na.rm=TRUE)
                  if (p == 0) { p <- 1 / (nperm + 1) } else { p <- p / nperm }
                },
                "less" = {
                  p <- sum(perms <= res, na.rm=TRUE)
                  if (p == 0) { p <- 1 / (nperm + 1) } else { p <- p / nperm }
                })
        res["p"] <- p
      }
      return(res)
    }

  }

  pvalsOrig = pvals
  if(is.null(pvals))
    pvals = sample(1:length(geneIds), length(geneIds))/length(geneIds)

  #convert geneEsts and pvals in case they are factors or character vectors
  geneEsts = as.numeric(as.character(geneEsts))
  pvals = as.numeric(as.character(pvals))

  geneDataDf = mapGenes(geneIds, geneEsts, pvals, forRankAndPlot = TRUE)

  if(is.null(pvalsOrig) & nrow(geneDataDf) != length(geneIds)){
    warning("repeat IDs were removed arbitrarily since p values were not supplied to base gene removal on.")
    pCut = FALSE
  }

  #deal with cmap rownames
  if(grepl("_at", rownames(drugPert)[1]))
    rownames(drugPert) = gsub("_at_", "", paste0(rownames(drugPert), "_"))
    #geneDataDf$ensemble = paste0(geneDataDf$ensemble, "_at")

  #shhould double check this against current L1000 dataset
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
  #changes here to use inputted featData if available/sent in
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

  if(cutOff > 1)
  {
    while(cutOff > 1)
      cutOff = cutOff/100

    warning("The cutOff value supplied was greater than 1 and has been reduced to a fraction between 0 and 1")
    if(pCut == TRUE & !is.null(pvalsOrig))
      print(paste("If it was intended that genes with a p value less than", cutOff, "should be kept, then the analysis will be conducted as intended"))
    #2 gene ests implies only direction is being used
    if(pCut == FALSE & length(unique(geneEsts)) > 2)
      print(paste("If it was intended that genes with an abs(estimate) in the top", cutOff*100, " percent should be kept, then the analysis will be conducted as intended"))
  }

  sigGenes = c(1:nrow(geneDataDf))
  if(pCut == TRUE & !is.null(pvalsOrig))
  {
    sigGenes = which(genesCmapAll$pvalue < cutOff)
    print(paste(length(sigGenes), "significant genes/genes with p <", cutOff,"in the analysis"))
  }else{
    if(length(unique(geneEsts)) > 2)
    {
      estOrd = order(abs(genesCmapAll$estimate), decreasing = TRUE)
      numKeep = floor(length(estOrd)*cutOff)
      sigGenes = estOrd[1:numKeep]
      print(paste(cutOff,"percent of genes with the top estimate will be used in the analysis for a total of", length(sigGenes), "genes"))
    }
  }

  genesCmap = genesCmapAll[sigGenes,]

  if(!is.na(extraData))
  {
    if(is.na(dim(extraData)))
    {
      extraData = cbind(extraData)
    }
    extraData = extraData[as.numeric(rownames(geneDataDf))]
    rowsToRem = removeGenes(idsOne = rownames(genesCmapAll), idsTwo = rownames(genesCmapAll), valuesTwo = extraData, conditionVals = extraCut, conditionsDirec = extraDirec)
    if(length(rowsToRem) > 0){
      genesCmap = genesCmap[-rowsToRem, ]
      print(paste("extraData was supplied and used to remove an additional ",length(rowsToRem) ," genes, the number of genes left in the analysis now is ", nrow(genesCmap)))
    }else
    {
      warning("no genes were left after using the criteria specified by extraData, extraCut, and extraDirec so the analysis will proceed without using the extra criteria")
    }
  }
  #when adding genes during analysis, add them in order of increasing p value. Don't know the nature of the estimates sent so safer to use p values
  #adding genes according to the order of estimates such as logFC or t-stat should not impact results much
  ordering = order(abs(genesCmap$pvalue), decreasing = FALSE)

  if(genesToStart >= 1)
  {
    while(genesToStart > 1)
      genesToStart = genesToStart/100

    warning("The genesToStart value supplied was greater than 1 and has been reduced to a fraction between 0 and 1")
    print(paste(genesToStart*100, "percent of the genes will be used in the first iteration of the analysis"))
    if(genesToStart == 1)
    {
      warning("As genesToSart was equal to 1 and would have crashed the program, it has been changed to 0.40")
      genesToStart = 0.40
    }
  }

  genesSigEnd = nrow(genesCmap)
  genesSigStart = ceiling(nrow(genesCmap)*genesToStart)
  genesSigInc = round((genesSigEnd - genesSigStart)/(numbIters - 1)) + 1

  print(paste0("Starting with ", genesSigStart, " genes and incrementing by ", genesSigInc, " genes. The last analysis will contain ", genesSigEnd," genes"))

  drugNotes = NULL
  if(!is.null(pharmSet))
    drugNotes = drugInfo(pharmSet)

  #########################################Section 2: Drug Ranking###########################################################

  quickTest = FALSE
  cmapList = list()
  naDrugRowsVec = c()
  pFrame = as.data.frame(NULL)
  fdrFrame = as.data.frame(NULL)
  csFrame = as.data.frame(NULL)
  print(paste("starting the",numbIters, "iterations"))
  for(g in 1:numbIters)
  {
    genesUse = genesSigStart + (g-1)*genesSigInc
    sigRows = ordering[1:genesUse]
    genesCmapUse = genesCmap[sigRows, ]

    cmapTab = as.data.frame(NULL)
    numDrugs = dim(drugPert)[2]
    drugNames = drugPert[1, ,1][1:numDrugs]
    drugNames = names(drugNames)
    noCores <- detectCores() - 2
    chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")

    if (nzchar(chk) && chk == "TRUE") {
      # use 2 cores in CRAN/Travis/AppVeyor
      noCores <- 2L
    } else {
      # use all cores in devtools::test()
      noCores <- parallel::detectCores()
    }
    noCores = 1
    # Initiate cluster
    cl <- makeCluster(noCores)
    pco = library(PharmacoGx)
    #clusterEvalQ(cl, library(PharmacoGx))

    if(quickTest == TRUE)
    {
      numDrugs = 300
      numbPerms = 100
    }

    drugsList = list()
    if(is.null(drugNameVec)) drugNameVec = colnames(drugPert)

    for(i in 1:length(drugNameVec))
    {
      drugColInd = which(colnames(drugPert) == drugNameVec[i])
      if(drugEst == TRUE)
      {
        drugsList[[i]] = drugPert[,drugColInd ,c("tstat", "pvalue")]
      }else{
        drugsList[[i]] = drugPert[,drugColInd ,c("estimate", "pvalue")]
      }
    }

    csFunc = function(z) connectivityScoreBroad(x = genesCmapUse, y = z,method=drugScoreMeth, gwc.method=gwcMethod, nperm = numbPerms)
    environment(csFunc) <- .GlobalEnv

    clusterExport(cl, list("drugsList", "genesCmapUse", "numbPerms", "drugScoreMeth", "gwcMethod", "connectivityScoreBroad", "gwc", "intersectList", "pco", "corWeighted", "combineTest"), envir = environment())
    dataCor = clusterApplyLB(cl, drugsList, csFunc)
    stopCluster(cl)

    for(i in 1:numDrugs)
    {
      if(!is.null(drugNotes))
      {
        drugInd = which(rownames(drugNotes) == drugNames[i])
        drugInfoTab = drugNotes[drugInd, ]
        cmapTab = rbind(cmapTab, c(dataCor[[i]], as.character(drugInfoTab)), stringsAsFactors = FALSE)
        colnames(cmapTab) = c("score", "pvalue", colnames(drugNotes[1,]))
      }else{
        cmapTab = rbind(cmapTab, c(dataCor[[i]]), stringsAsFactors = FALSE)
        colnames(cmapTab) = c("score", "pvalue")
      }
    }

    rownames(cmapTab) = drugNames[1:numDrugs]
    #sorting = order(cmapTab[,"score"], decreasing = TRUE)
    #cmapTabSort = cmapTab[sorting,]
    cmapTabNew = cbind(cmapTab[, 1:2], p.adjust(cmapTab[,2], method = "fdr"), cmapTab[, 3:ncol(cmapTab)])
    colnames(cmapTabNew)[3] = "fdr adjusted p"
    cmapTab = cmapTabNew


    #one of the drugs score kept showing up as NA, so remove it
    naDrugRows = which(is.na(cmapTab[,1]))
    naDrugRowsVec = c(naDrugRowsVec, naDrugRows)
    #if(length(naDrugRows) > 0)
    #  cmapTab = cmapTab[-naDrugRows, ]

    #rbind(paste(as.character(vec), collapse = ", "))
    if(g == 1)
    {
      pFrame = cbind(format(as.numeric(as.character(cmapTab[, "pvalue"])), scientific =  TRUE, digits = 3))
      fdrFrame = cbind(format(as.numeric(as.character(cmapTab[, "fdr adjusted p"])), scientific =  TRUE, digits = 3))
      csFrame = cbind(sprintf("%.5f", as.numeric(as.character(cmapTab[, "score"]))))
    }else{
      pFrame = cbind(pFrame, format(as.numeric(as.character(cmapTab[, "pvalue"])), scientific =  TRUE, digits = 3))
      fdrFrame = cbind(fdrFrame, format(as.numeric(as.character(cmapTab[, "fdr adjusted p"])), scientific =  TRUE, digits = 3))
      csFrame = cbind(csFrame, sprintf("%.5f", as.numeric(as.character(cmapTab[, "score"]))))
    }


    cmapList[[g]] = cmapTab
    #write.table(cmapTabSort, paste("genes used is", genesUse, fileNameCmap), col.names = NA, sep="\t")
    print(g)
  }
  rownames(pFrame) = drugNames[1:numDrugs]
  rownames(fdrFrame) = drugNames[1:numDrugs]
  rownames(csFrame) = drugNames[1:numDrugs]

  if(length(naDrugRowsVec) > 0)
  {
    naDrugRowsVec = unique(naDrugRowsVec)
    pFrame = pFrame[-naDrugRowsVec, ]
    csFrame = csFrame[-naDrugRowsVec, ]
    fdrFrame = fdrFrame[-naDrugRowsVec, ]
    for(g in 1:numbIters)
    {
      cmapList[[g]] = cmapList[[g]][-naDrugRowsVec, ]
    }
  }
  #rank the drugs according to their connectivity score and compute the change in their rank over the iterations
  rankFrame = as.data.frame(NULL)
  #below frame finds drugs that mimick the phenotype the most
  rankFrameRev = as.data.frame(NULL)
  rankChangeFrame = as.data.frame(NULL)
  colnameVec = c()
  for(g in 1:numbIters)
  {
    ranking = c()
    rankingRev = c()
    cmapRes = cmapList[[g]]["score"]
    ordering = sort(-1*as.numeric(cmapRes[,1]), index.return=TRUE,decreasing = TRUE)$ix
    orderingRev = sort(as.numeric(cmapRes[,1]), index.return=TRUE,decreasing = TRUE)$ix
    for(i in 1:length(ordering))
    {
      ranking[ordering[i]] = i
      rankingRev[orderingRev[i]] = i
    }
    if(g == 1)
    {
      rankFrame = as.data.frame(ranking)
      rankFrameRev = as.data.frame(rankingRev)
    }
    if(g > 1)
    {
      rankFrame = cbind(rankFrame, ranking)
      rankFrameRev = cbind(rankFrameRev, rankingRev)
    }
    if(g == 2)
      rankChangeFrame = as.data.frame(abs(ranking - oldRanking))
    if(g > 2)
      rankChangeFrame = cbind(rankChangeFrame, abs(ranking - oldRanking))
    colnameVec = c(colnameVec, as.character(genesSigStart + (g-1)*genesSigInc))
    oldRanking = ranking
  }
  rownames(rankFrame) = rownames(cmapRes)
  colnames(rankFrame) = colnameVec

  rownames(rankFrameRev) = rownames(cmapRes)
  colnames(rankFrameRev) = colnameVec

  rownames(rankChangeFrame) = rownames(cmapRes)
  colnames(rankChangeFrame) = colnameVec[2:length(colnameVec)]

  rankChangeFrameOrig = rankChangeFrame
  for(g in 1:dim(rankChangeFrame)[2])
  {
    changes = rankChangeFrame[,g]
    changeRanks = rank(-1*changes, ties.method = "min")
    changeRanksOrig = rank(changes, ties.method = "min")
    rankChangeFrame[,g] = changeRanks
    rankChangeFrameOrig[,g] = changeRanksOrig
  }

  #if a drug ranked higher on the list that ranked drugs according to negative connectivity scores use the score from that list as its final
  #score and multiply the score by a negative 1. Otherwise use the rank from other (positive connectivity score) list
  drugRankMeansPos = rowMeans(rankFrame)
  drugRankMeansNeg = rowMeans(rankFrameRev)
  drugRankSigns = c()
  drugRankMeans = c()
  for(g in 1:length(drugRankMeansPos))
  {
    if(drugRankMeansPos[g] >= drugRankMeansNeg[g])
    {
      drugRankMeans[g] = drugRankMeansPos[g]
      drugRankSigns[g] = 1
    }
    if(drugRankMeansNeg[g] > drugRankMeansPos[g])
    {
      drugRankMeans[g] = drugRankMeansNeg[g]
      drugRankSigns[g] = -1
    }
  }

  drugRankChangeMeans = rowMeans(rankChangeFrame)
  drugRankChangeMeansOrig = rowMeans(rankChangeFrameOrig)
  #final score is an average of the drugs mean rank and the rank of its stability from iteration to iteration
  sigFinalScores = drugRankSigns*(drugRankMeans + drugRankChangeMeans)/2

  pCol = c()
  fdrCol = c()
  csCol = c()
  for(i in 1:nrow(csFrame))
  {
    pCol[i] = paste(pFrame[i, ], collapse = ", ")
    fdrCol[i] = paste(fdrFrame[i, ], collapse = ", ")
    csCol[i] = paste(csFrame[i, ], collapse = ", ")
  }
  #although best final scores for reversing the phenotypes are large negative ranks, paste drugRankMeansPos and drugRankMeansOrig to the table
  #as they are easiest to interpret. A value of 1 on these lists means the drug had the largest negative connectivity score each iteration (rankMeans)
  #and its rank changed the least from iteration to iteration relative to all the other drugs (rankChange)
  cmapResults = cbind(sigFinalScores, drugRankMeansPos, drugRankChangeMeansOrig, csCol, pCol, fdrCol, cmapList[[1]][4:ncol(cmapList[[1]])])
  sorting = order(cmapResults[,"sigFinalScores"], decreasing = FALSE)
  cmapResults = cmapResults[sorting,]
  colnames(cmapResults)[1:6] = c("Final Score", "Mean Rank (1 implies most negative connectivity scores)", "Rank of Mean Rank Change (lower implies more stable)", "connectivity scores", "p values", "fdr adjusted p")

  #rankFrameInd = which(rownames(rankFrame) == topDrug)
  #rankChangeInd = which(rownames(rankChangeFrame) == topDrug)
  #rankVsSize = rankFrame[rankFrameInd, ]
  #changeVsSize = rankChangeFrameOrig[rankChangeInd,]

  #if(sum(grepl("BRD-", rownames(cmapResults))) > nrow(cmapResults)/2)
  #{
  #  rownames(cmapResults) = paste(cmapResults$pert_iname, " (ID ", rownames(cmapResults),")", sep = "")
  #  #give drugPert appropr
  #  for(i in 1:length(colnames(drugPert)))
  #    colnames(drugPert)[i] = paste(drugNotes$pert_iname[which(colnames(drugPert)[i] == rownames(drugNotes))], " (ID ",colnames(drugPert)[i],")", sep = "")

  #}
  #topDrug = rownames(cmapResults)[1]

#########################################Section 3: Results Validation###########################################################

  #gsea ranks vs genes seem to trend upwards, add genes previous best ranked drug has higher ranks, gwc seems to stabilize
  if(inspectResults == TRUE)
  {
    #sizeVec = c()
    #for(g in 1:numbIters)
    #  sizeVec = c(sizeVec, genesSigStart + (g-1)*genesSigInc)

    #plot(sizeVec, rankVsSize, main = paste("Rank as a Function of Query Size for the Top Drug", topDrug), xlab = "Query Size (number of genes in the analysis)", ylab = "Rank (1 implies best at reversing disease phenotype)")
    #plot(1:(numbIters-1), changeVsSize, main = paste("Change in Rank from Query to Query for the Top Drug", topDrug), xlab = "Query Number", ylab = "abs(Rank of Mean Rank Change) Between Queries")

    direcOnly = FALSE
    if(is.null(pvalsOrig) | length(unique(geneEsts)) == 2)
      direcOnly = TRUE

    topDrug = rownames(cmapResults)[1]
    #topDrug = drugName
    rankFrameInd = which(rownames(rankFrame) == topDrug)
    rankChangeInd = which(rownames(rankChangeFrame) == topDrug)
    #should probably mae sure these rank frames and ran change frame are returned, or the results are in returned frame
    rankVsSize = rankFrame[rankFrameInd, ]
    changeVsSize = rankChangeFrameOrig[rankChangeInd,]

    dev.new()
    plot(1:length(rankVsSize), rankVsSize, main = paste("Rank as a Function of the Iterations for", topDrug), xlab = "Query Number", ylab = "Rank (1 implies best at reversing disease phenotype)")
    dev.new()
    plot(1:length(changeVsSize), changeVsSize, main = paste("Change in Rank from Query to Query for", topDrug), xlab = "Query Number", ylab = "abs(Rank of Mean Rank Change) Between Queries")

    #drugRankPlots(cmapResults, topDrug)

    if(is.na(volcPlotEsts[1]) | is.null(pvalsOrig)){
      volcPlotEsts = geneDataDf$geneEsts
    }else{
      volcPlotEsts = volcPlotEsts[as.numeric(rownames(geneDataDf))]
    }

    idsUsed = rownames(genesCmap)
    idsUsed = which(geneDataDf$ensemble %in% idsUsed)
    #remove inspectDrugResults and just put code below here

    genesEst = volcPlotEsts[idsUsed]
    genesSymb = geneDataDf$symbol[idsUsed]
    genesP = geneDataDf$pvals[idsUsed]
    genesId = geneDataDf$ensemble[idsUsed]
    drugResults = cmapResults
    topDrugName = rownames(drugResults)[1]


    if(direcOnly == FALSE){
      #dev.new()
      makeVolcPlot(genesEst, genesSymb, genesP, genesId, rownames(drugPert))
    }

    #below not currently working for fgsea and gwcCmapBox
    driveGenesNegScore = plotDrugGeneInteraction(genesEst, genesSymb, genesP, genesId, drugPert, topDrugName, drugScoreMeth, gwcMethod, droveNegScore = TRUE, drugEst, drugVolcPlot = FALSE, supressPlot = FALSE)
    if(drugScoreMeth == "gwc" | drugScoreMeth == "fgsea")
      print(paste("For the plots that show the genes that drove the negative connectivity score, blue and red labelled genes should cross 0 when going from the data plot to drug", topDrugName,"volcano plot and the blue dots in the plot for the data should be genes that have lower expression in the disease state relative to the normal state. If the prior statements are not true then these are not the drugs you are looking for. You likely need to switch the sign of the supplied geneEsts vector."))

    genesOfInt = c(driveGenesNegScore$experimentList$posGeneDriversExper$geneSymb, driveGenesNegScore$experimentList$negGeneDriversExper$geneSymb)
    geneChangesFromDrug = getDrugsImpactOnGenes(drugPert, pharmSet, mDataType, topDrugName, genesOfInt)
    if(drugScoreMeth == "gwc" | drugScoreMeth == "fgsea")
      print(paste("If the negative gene estimate genes from the data (blue dots in data plots) that drove the negative connectivity score have not had their expression levels increase (Percent Change > 0 and usually blue boxes) on the heatmaps for drug", topDrugName, "then this drug does not appear to increase the expression of genes that had low expression in the disease phenotype relative to the normal phenotype and this is not the drug you are looking for (blue dots in plots should be blue in the heatmap for the drug that reverses the phenotype). "))

    if(showMimic == TRUE)
    {
      botDrugName = rownames(drugResults)[nrow(drugResults)]
      driveGenesPosScore = plotDrugGeneInteraction(genesEst, genesSymb, genesP, genesId, drugPert, botDrugName, drugScoreMeth, gwcMethod, droveNegScore = FALSE, drugEst, drugVolcPlot = FALSE, supressPlot = FALSE)
      genesOfIntBot = c(driveGenesPosScore$experimentList$posGeneDriversExper$geneSymb, driveGenesPosScore$experimentList$negGeneDriversExper$geneSymb)

      geneChangesFromBotDrug = getDrugsImpactOnGenes(drugPert, pharmSet, mDataType, botDrugName, genesOfIntBot)
    }

  }

  #add clinical trial info if available
  data("cmapTrialInfo")
  cmapNames = tolower(rownames(cmapResults))
  repoNames = tolower(cmapTrialInfo$drug_name)

  if(length(which(cmapNames %in% repoNames)) > 0)
  {
    trialInfoVec = c()
    trialStatusVec = c()
    for(i in 1:nrow(cmapResults))
    {
      trialInds = which(tolower(rownames(cmapResults)[i]) == tolower(cmapTrialInfo$drug_name))
      if(length(trialInds) > 0){
        drugClinInfo = cmapTrialInfo[trialInds, ]
        trialInfo = drugClinInfo$Indication
        trialInfoVec[i] = paste(trialInfo, collapse = " | ")
        trialStatus = drugClinInfo$status
        trialStatusVec[i] = paste(trialStatus, collapse = " | ")
      }else{
        trialInfoVec[i] = NA
        trialStatusVec[i] = NA
      }
    }
    cmapResults = cbind(cmapResults, trialInfoVec, trialStatusVec)
    colnames(cmapResults)[(ncol(cmapResults) - 1):(ncol(cmapResults))] = c("Drug Clinical Trials", "Trial Outcomes")
  }
  #double check that each method returns connectivity scores between 0 and 1!

  return(cmapResults)
}
#to add new drug repurposing technique
#1. add technique to connectivityScore function and ensure more negative outcome implies reverses the phenotype (for ranking later)
#2. adjust getSigGeneVolcPlot to be able to determine which genes drove the connectivity score using that technique
#3. add statement to inspectDrugResults function explaining how to interpret data to determine if the analysis is correct
#4. adjust drugScoreMeth in file info section to include the technique
#5. adjust vignette or other areas that should make reference to the technique

#1250 genes at iteration 10 in original data received

#library("devtools")
#library(roxygen2)
#setwd("C:\\Users\\micha\\Documents\\PMH Research\\gwcCMapCur")
#document()
#setwd("..")
#install("gwcCMapCur")
#library(gwcCMapCur)

#for testing the functions
#data("geneDataGwc")
#data("drugPertEx")
#data("psetSub")
#geneIds = geneDataGwc$symbol
#geneEsts = geneDataGwc$t
#pvals = geneDataGwc$P.Value
#drugPert = drugPertEx
#pharmSet = psetSub
#drugScoreMeth = "gwc"
#volcPlotEsts = geneDataGwc$logFC
#pCut = TRUE
#cutOff = 0.05
#genesToStart = 0.20
#numbIters = 10
#gwcMethod = "spearman"
#numbPerms = 1000
#drugEst = TRUE
#inspectResults = TRUE
#showMimic = FALSE
#mDataType = "rna"
#extraData = NA
#extraCut = NA
#extraDirec = NA
#drugNameVec = NULL
