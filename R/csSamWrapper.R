#' csSamWrapper function - performs entire functionality
#' 
#' csSamWrapper function - performs entire functionality
#' 
#' @param G Matrix of gene expression, columns ordered in the same order at the
#' cell-frequency matrix (n by g, n samples, g genes)
#' @param cc Matrix of cell-frequency. (n by k, n samples, k cell-types)
#' @param y A numeric vector of group association of each sample. Either 1 or 2.
#' @param nperms The number of permutations to perform.
#' @param alternative two.sided less greater
#' @param standardize Standardize sample or not. Default is TRUE.
#' @param medianCenter Median center rhat distributions. Default is TRUE.
#' @param logRm Exponentiate data for deconvolution stage. Default is FALSE
#' @param logBase Base of logaritm used to determine exponentiation factor.
#' Default is 2
#' @param nonNeg For single channel arrays. Set any cell-specific expression
#' estimated as negative, to a ceiling of 0. It is conservative in its study of
#' differential expression. Default is TRUE.
#' @param fileName PDF file containing plots of FDR vs. number of genes called
#' for whole tissue comparison (via SAM) as well as each cell-type (by csSAM)
#' @return Returns a list containing:
#' \item{csSAM}{A list object containing a fit (cell-type specfic
#' expression) for each group. Each element in the list is an object returned by
#' csFit.}
#' \item{SAM}{A list output of the fdrSAM function.}
#' \item{sigGene.csSAM}{A list of significant genes.}
#' \item{fileName}{The filename into whcih the FDR plots are dumped.}
#' @author Shai Shen-Orr, Rob Tibshirani, Narasimhan Balasubramanian, David Wang
#' @seealso
#' \code{\link{csfit}},\code{\link{csSAM}},\code{\link{fdrCsSAM0}},\code{\link{plotCsSAM}}
#' @cite Shen-Orr2010
#' @export 
#' 
csSamWrapper <- function(G,cc,y,nperms = 200,alternative = 'two.sided',standardize=TRUE,medianCenter=TRUE, logRm =FALSE,logBase = 2,nonNeg=TRUE,fileName='csSAMout.pdf') {

  # run csSAM
  csfit <- csSAMfit(x = t(G), cc = t(cc), y = y, nperms = nperms
            , alternative = alternative, standardize = standardize, medianCenter = medianCenter
            , logRm = logRm, logBase = logBase, nonNeg = nonNeg)
	# run SAM
  tt.sam <- runSAM(G, y)
	fdr.sam <- fdrSAM(G, y, nperms=nperms, tt.sam, alternative)
  sigGene <- csTopTable(csfit, alternative = alternative, merge = TRUE)

	plotCsSAM(csfit, fdr.sam, alternative, fileName = fileName)
	return(list(csSAM = basisfit(csfit), SAM = fdr.sam, sigGene.csSAM = sigGene, fileName = fileName))
  
}
