#'Function to calculate a modified gwc score between two vectors
#'
#'As opposed to calculating a genome wide correlation as is done by PharmacoGx::gwc, this function calculates a modifed gwc
#'by only using the sign of the second vector and not weighting the p values of the second vector. The rationale behind this is 
#'that other drug repurposing techniques only use the signs of the gene signatures (ex gsea) supplied in order to identify drug candidates.
#'
#'@examples
#'data(CCLEsmall)
#'x <- molecularProfiles(CCLEsmall,"rna")[,1]
#'y <- molecularProfiles(CCLEsmall,"rna")[,2]
#'x_p <- rep(0.05, times=length(x))
#'y_p <- rep(0.05, times=length(y))
#'names(x_p) <- names(x)
#'names(y_p) <- names(y)
#'gwc(x,x_p,y,y_p, nperm=100)
#'
#'@param x1 \code{numeric} vector of effect sizes (e.g., fold change or t statitsics) for the first experiment
#'@param p1 \code{numeric} vector of p-values for each corresponding effect size for the first experiment
#'@param x2 \code{numeric} sig of the effect size (e.g., fold change or t statitsics) for the second experiment
#'@param method.cor \code{character} string identifying if a \code{pearson} or
#'\code{spearman} correlation should be used
#'@param nperm \code{numeric} how many permutations should be done to determine the p value
#'@param truncate.p \code{numeric} Truncation value for extremely low p-values
#'
#'@return \code{numeric} a vector of two values, the correlation and associated p-value.
#'@export
##            -
##
## http://en.wikipedia.org/wiki/Pearson_product-moment_correlation_coefficient#Calculating_a_weighted_correlation
## http://www.mathworks.com/matlabcentral/fileexchange/20846-weighted-correlation-matrix
## F. Pozzi, T. Di Matteo, T. Aste, "Exponential smoothing weighted correlations", The European Physical Journal B, Vol. 85, No 6, 2012. DOI: 10.1140/epjb/e2012-20697-x
#################################################


gwcCmapBox = function (x1, p1, x2, method.cor = c("pearson", "spearman"), nperm = 10000, truncate.p = 1e-16) 
{
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