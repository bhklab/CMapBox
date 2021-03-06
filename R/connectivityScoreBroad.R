

########################
## Benjamin Haibe-Kains & Cat Olsen
## September 9, 2014
########################



#' Function computing connectivity scores between two signatures
#'
#' A function for finding the connectivity between two signatures, using either
#' the GSEA method based on the KS statistic, or the gwc method based on a
#' weighted spearman statistic. The GSEA analysis is implemented in the piano package.
#'
#' @references
#'    F. Pozzi, T. Di Matteo, T. Aste, "Exponential smoothing weighted
#'    correlations", The European Physical Journal B, Vol. 85, No 6, 2012. DOI:
#'    10.1140/epjb/e2012-20697-x
#' @references
#'    Varemo, L., Nielsen, J. and Nookaew, I. (2013) Enriching the gene set
#'    analysis of genome-wide data by incorporating directionality of gene
#'    expression and combining statistical hypotheses and methods. Nucleic
#'    Acids Research. 41 (8), 4378-4391. doi: 10.1093/nar/gkt111
#' @references
#'    Cheng, Jie, et al. "Systematic evaluation of connectivity map for disease
#'    indications." Genome medicine 6.12 (2014): 95.
#'
#' @examples
#' xValue <- c(1,5,23,4,8,9,2,19,11,12,13)
#' xSig <- c(0.01, 0.001, .97, 0.01,0.01,0.28,0.7,0.01,0.01,0.01,0.01)
#' yValue <- c(1,5,10,4,8,19,22,19,11,12,13)
#' ySig <- c(0.01, 0.001, .97,0.01, 0.01,0.78,0.9,0.01,0.01,0.01,0.01)
#' xx <- cbind(xValue, xSig)
#' yy <- cbind(yValue, ySig)
#' rownames(xx) <- rownames(yy) <- c('1','2','3','4','5','6','7','8','9','10','11')
#' data.cor <- connectivityScore(xx,yy,method="gwc", gwc.method="spearman", nperm=300)
#'
#' @param x A \code{matrix} with the first gene signature. column one being the values for the genes
#' and column two being their corresponding significance values.
#' @param y A \code{matrix} with the second signature. column one being the values for the genes
#' and column two being their corresponding significance values.
#' @param method \code{character} string identifying which method to use, out of 'gsea', 'xsum' and 'gwc'
#' @param nperm \code{numeric}, how many permutations should be done to determine
#'   significance through permutation testing? The minimum is 100, default is
#'   1e4.
#' @param nthread \code{numeric}, how many cores to run parallel processing on.
#' @param gwc.method \code{character}, should gwc use a weighted spearman or pearson
#'   statistic?
#' @param ... Additional arguments passed down to gsea and gwc functions
#' @return \code{numeric} a numeric vector with the score and the p-value associated
#'   with it
#' @export
#'
#'

connectivityScoreBroad <- function(x, y, method=c("fgsea", "gwc", "xsum", "gwcCmapBox"), nperm=1e4, nthread=1, gwc.method=c("spearman", "pearson"), ...) {

  method <- match.arg(method)
  if (class(x) != "matrix") {
    x <- as.matrix(x)
  }
  if (class(y) != "matrix") {
    y <- as.matrix(y)
  }
  if ((ncol(x) != 2 || ncol(y) != 2) && method=="gwc") {
    stop ("x and y should have 2 columns: effect size and corresponding p-values")
  }

  if (method == "fgsea" && nrow(y) >= nrow(x)) {
    warning("GSEA method: query gene set (y) larger than signature (x)")
  }

  if (is.null(rownames(x)) || is.null(rownames(y)) || !length(intersect(rownames(x), rownames(y)))) {
    stop ("Row names of x and y are either missing or have no intersection")
  }
  if (nperm < 100){
    stop ("The minimum number of permutations for permutation testing is 100")
  }
  switch (method,
          "fgsea" = {
            ## remove missing values
            y <- y[!is.na(y[ ,1]), , drop=FALSE]
            x <- x[!is.na(x[ ,1]), , drop=FALSE]
            ## create gene set
            gset <- cbind("gene"=rownames(y), "set"=ifelse(as.numeric(y[ , 1]) >= 0, "UP", "DOWN"))
            gset <- piano::loadGSC(gset)
            ## run enrichment analysis
            nes <- piano::runGSA(geneLevelStats=x[ , 1], geneSetStat="fgsea", gsc=gset, nPerm=nperm + (nperm %% nthread), ncpus=nthread, verbose=FALSE, adjMethod="none",...)
            ## merge p-values for negative and positive enrichment scores
            nes$pDistinctDir <- nes$pDistinctDirUp
            nes$pDistinctDir[is.na(nes$pDistinctDirUp), 1] <- nes$pDistinctDirDn[is.na(nes$pDistinctDirUp), 1]
            nes.up <- c(nes$statDistinctDir[which(names(nes$gsc) == "UP"), 1], nes$pDistinctDir[which(names(nes$gsc) == "UP"), 1])
            nes.down <- c(nes$statDistinctDir[which(names(nes$gsc) == "DOWN"), 1], nes$pDistinctDir[which(names(nes$gsc) == "DOWN"), 1])
            ## combine UP and DOWN
            #divide vy 2?
            if (length(nes.up) == 0){
              score = c("es" = -nes.down[1], "p" = nes.down[2])
            } else if (length(nes.down) == 0){
              score = c("es" = nes.up[1], "p" = nes.up[2])
            } else if (complete.cases(cbind(nes.up, nes.down)) && sign(nes.up[1]) != sign(nes.down[1])) {
              score <- c("es"=(nes.up[1] - nes.down[1]) / 2, "p"=combineTest(p=c(nes.up[2], nes.down[2]), method="fisher", na.rm=TRUE))
            } else {
              score <- c("score"=0, "p"=1)
            }
          },
          "gwc" = {
            ## intersection between x and y
            ii <- intersect(rownames(x), rownames(y))
            if(length(ii) < 10) {
              stop ("Less than 10 probes/genes in common between x and y")
            }
            score <- gwc(x1=x[ii, 1], p1=x[ii, 2], x2=y[ii, 1], p2=y[ii, 2], method.cor=gwc.method, nperm=nperm, ...)
            names(score) <- c("score", "p")
          },
          "gwcCmapBox" = {

            gwcCmapBox = function(x1, p1, x2, method.cor = c("pearson", "spearman"), nperm = 10000, truncate.p = 1e-16)
            {
              #corWeighted below shouldn't be needed and can probably be imported from pharmacoGx
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

              method.cor = match.arg(method.cor)
              ii = intersectList(names(x1), names(p1), names(x2))
              if (length(ii) < 10) {
                stop("Less than 10 probes/genes in common between x and y")
              }
              x1 = x1[ii]
              p1 = p1[ii]
              x2 = sign(x2[ii])
              p1[!is.na(p1) & p1 < truncate.p] <- truncate.p
              p1 = -log10(p1)
              p1 = p1/sum(p1, na.rm = TRUE)
              w = p1
              res = corWeighted(x = x1, y = x2, w = w, method = method.cor,
                                nperm = nperm)
              return(res)
            }
            ## intersection between x and y
            ii <- intersect(rownames(x), rownames(y))
            if(length(ii) < 10) {
              stop ("Less than 10 probes/genes in common between x and y")
            }
            score <- gwcCmapBox(x1=y[ii, 1], p1=y[ii, 2], x2=x[ii, 1], method.cor=gwc.method, nperm=nperm)
            names(score) <- c("score", "p")

          },
          "xsum" = {
            ii <- intersect(rownames(x), rownames(y))
            if(length(ii) < 10) {
              stop ("Less than 10 probes/genes in common between x and y")
            }
            #look over with haibe-kains
            #Q trust magnitudes in drugPert or rank?
            #Q 1- p if very positive score?
            #Method is different as it gives best -ve score to compounds that have significant genes in disease as the most
            #effected genes by the compound (effect in favourable direction)
            #paper takes top N up in disease, top N down in disease, top N up in compound, top N down in compound
            #in order to use all significant genes, slight adjustments to papers method
            #1. use all genes sent in for disease, so top N1 up in disease and top N2 down in disease with numGenes = N1 + N2
            #2. implies look for those top N1 genes in compound and sum them, then seek the bot N2 genes in compound and sum them
            #   score = sum(N1 genes present) - sum(N2 genes present) as one wants to reverse the phenotype
            #3.. create permutation test to get p value
            y1 = y[, 1]
            names(y1) = rownames(y)
            x1 = x[, 1]
            names(x1) = rownames(x)
            sumUpX = 0
            sumDownX = 0
            upInDis = x1[which(x1 >= 0)]
            downInDis = x1[which(x1 < 0)]
            topAndBotCompound = c(y1[order(y1)][1:length(upInDis)], y1[order(y1, decreasing = TRUE)][1:length(downInDis)])
            upInDisInter = intersect(names(upInDis), names(topAndBotCompound))
            downInDisInter = intersect(names(downInDis), names(topAndBotCompound))

            if(length(upInDisInter) > 0)
              sumUpX = sum(topAndBotCompound[upInDisInter])
            if(length(downInDisInter) > 0)
              sumDownX = sum(topAndBotCompound[downInDisInter])
            scoreVal = sumUpX - sumDownX

            randScoreVec = c()
            if((length(upInDisInter) + length(downInDisInter)) > 0){
              for(i in 1:nperm)
              {
                sumUpXRand = 0
                sumDownXRand = 0
                #possible occasional gene overlap in sumUpXRand and sumDownXRand but impact should be inconsequential
                upInDisRand = sample(1:length(topAndBotCompound), length(upInDisInter))
                if(length(upInDisRand) > 0)
                  sumUpXRand = sum(topAndBotCompound[upInDisRand])
                downInDisRand = sample(1:length(topAndBotCompound), length(downInDisInter))
                if(length(downInDisRand) > 0)
                  sumDownXRand = sum(topAndBotCompound[downInDisRand])
                randScore = sumUpXRand - sumDownXRand
                randScoreVec = c(randScoreVec, randScore)
              }
              p = length(which(randScoreVec <= scoreVal))/nperm
              if(p == 0)
                p = 1/nperm
            }else{
              p = 1
            }
            score = c(scoreVal, p)
            names(score) <- c("score", "p")
          }
  )
  return (score)
}

