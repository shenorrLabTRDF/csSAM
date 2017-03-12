#' findSigGene
#' 
#' Find the false discovery rate for each gene in each cell-type.
#' 
#' 
#' @note From version 1.3, the FDR computation is implemented in \emph{C++}, 
#' using \pkg{Rcpp}.
#'
#' @param rhat Matrix of cell-specific contrasts for each gene in each cell-type
#' as computed for the original group classification.
#' @param csSAMData List object returned from fdrCsSAM.
#' @param alternative alternative hypothesis being tested: \code{'two.sided'}, \code{'less'} 
#' or \code{'greater'}
#' @param version number that indicates the version of computation to perform: 1 or 2.
#' Current default is \code{version = 2}, which fixes an issue for alternatives 
#' \code{'less'} and \code{'greater'}. Use \code{version = 1} to obtain old behaviour.
#' @return A matrix size k by g where k is the number of cell-types and g is the
#' number of genes. For each cell in the matirx, listed is the FDR of the gene
#' for a difference in a given cell-type.
#' @author Shai Shen-Orr, Rob Tibshirani, Narasimhan Balasubramanian, David Wang
#' 
#' \emph{C++} implementation/optimisation by Renaud Gaujoux 
#' @cite Shen-Orr2010
#' @useDynLib csSAM
findSigGene <- function(rhat, csSAMData, alternative = c('two.sided', 'less', 'greater'), version = getOption('csSAM.siggene.version', 2)) {
	
	# version 1 only computed two.sided
	if( version == 1 ) alternative <- 'two.sided'
	alt <- match.arg(alternative)
	
	if( alt == 'less' ) rhat <- -rhat
	else if( alt == 'two.sided' ) rhat <- abs(rhat)
	
	.Call('findSigGenes', rhat, csSAMData$cutp.g, csSAMData$fdr.g, PACKAGE='csSAM')
}


## old plain R version
.findSigGene_R <- function(rhat, csSAMData, alternative = c('two.sided', 'less', 'greater'), version = getOption('csSAM.siggene.version', 2)) {
	#numgene=ncol(G)
	numgene = ncol(rhat)
	#numcell = ncol(cc)
	numcell = nrow(rhat)
	thresholdVec = csSAMData$fdr.g
	cutoff <- csSAMData$cutp.g
	thresholdLen = length(thresholdVec[numcell,])
	sigGene <- array(dim = c(numcell, numgene))
	sigGene[,] = 1
	
	# version 1 only computed two.sided
	if( version == 1 ) alternative <- 'two.sided'
	alt <- match.arg(alternative)
	
	if( alt == 'less' ) rhat <- -rhat
	else if( alt == 'two.sided' ) rhat <- abs(rhat)
	
	for (curThresh in 1:thresholdLen) {
		for (curcell in 1:numcell) {
			for (curgene in 1:numgene) {
				if( rhat[curcell,curgene] >= cutoff[curcell,curThresh]) {
					sigGene[curcell,curgene] = thresholdVec[curcell,curThresh]
				}
			}
		}
	}
	
	return (sigGene)
}

#.findSigGene_cmp <- compiler::cmpfun(.findSigGene_R) 
