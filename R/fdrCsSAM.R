#' fdrCsSAM - Deprecated Version
#' 
#' Estimates the false discovery rate for the identified cell-specific
#' differences in gene expression.
#' 
#' @note From version 1.3, the FDR computation is implemented in \emph{C++}, 
#' using \pkg{Rcpp}.
#' 
#' @param G Matrix of gene expression, columns ordered in the same order at the
#' cell-frequency matrix (n by p, n samples, p genes)
#' @param cc Matrix of cell-frequency. (n by k, n samples, k cell-types)
#' @param y A numeric vector of group association of each sample. Either 1 or 2.
#' @param n A nuermic vector describing the number of samples in a group
#' @param numcell The number of cell-types to consider
#' @param numgene The number of genes being considered
#' @param rhat The contrast in cell-type expression for each cell-type as
#' observed between the two groups being compared.
#' @param nperms The number of permutations to perform.
#' @param alternative Type of test to conduct - choose between
#' 'two.sided','greater',or 'less'
#' @param standardize Standardize sample or not. Default is TRUE
#' @param medianCenter Median center rhat distributions. Default is TRUE.
#' @param logRm Exponentiate data for deconvolution stage. Default is FALSE
#' @param logBase Base of logaritm used to determine exponentiation factor.
#' Default is 2
#' @param nonNeg For single channel arrays. Set any cell-specific expression
#' estimated as negative, to a ceiling of 0. It is conservative in its study of
#' differential expression. Default is FALSE.
#' @return A list.
#' \item{fdr.g}{A matirx false dicovery rates for csSAM comparison for each
#' cell-type at different thresholds. A set of 100 theresholds is determined
#' automatically from the data (k by 100, where k is number of cells).}
#' \item{avrhatperm}{}
#' \item{rhatperm}{A matrix sized pXkXg which stores the contrast of a
#' given gene g in cell type k in permutation p of the data.}
#' \item{cutp.g}{A matrix k by 100, where k is the number of cell tpes.
#' Lists the 100 cutoff thresholds for each cell-type as determined
#' automatically from the computed contrast.}
#' \item{rhat}{A matrix object with the result of contrasting the average
#' cell-specific expression profile of the two groups, per cell-type (Size k by
#' g where k is the number of cells and g is the number of genes).}
#' \item{ncall.g}{Number of genes called significant at the given cutoff
#' threshold with a FDR matching that indicated in fdr.g}
#' @author Shai Shen-Orr, Rob Tibshirani, Narasimhan Balasubramanian, David Wang
#' 
#' \emph{C++} implementation/optimisation by Renaud Gaujoux.
#' @cite Shen-Orr2010
fdrCsSAM0 <-
function (G,cc,y,n,numcell,numgene,rhat,nperms,alternative=c('two.sided', 'greater', 'less')
        ,standardize=TRUE,medianCenter=TRUE,logRm=FALSE,logBase = 2,nonNeg=FALSE) {

  .Deprecated('fdrCsSAM')
    alternative <- match.arg(alternative)
  numgene=ncol(G)
  rhatperm <- array(dim = c(nperms,numcell,numgene))
  perm = list()
  
  for (i in 1:nperms) {
    o=sample(1:length(y))
    ystar = y[o]
    for (curset in 1:2) {
      perm[[curset]]= csfit(cc[ystar==curset,], G[ystar==curset,],logRm,logBase)
    }
    rhatperm[i,,] = csSAM(perm[[1]]$ghat, perm[[1]]$se, n[1], perm[[2]]$ghat, perm[[2]]$se, n[2],
              standardize, medianCenter, nonNeg)
  }
#  cutp.g=matrix(NA,nrow=numcell,ncol=100)
#  numcut = ncol(cutp.g)
#  
#  fdr.g=ncall.g=nperm.g<-array(dim = c(numcell, numcut))
#  
#  for(j in 1:numcell)
#    cutp.g[j,]=seq(0,max(abs(rhat[j,])),length=100)	
#  for (i in 1:numcut) {
#    for (curcell in 1:numcell) {
#		if(alternative == 'two.sided') {
#			fdr.g[curcell,i]=sum(abs(rhatperm[,curcell,])>cutp.g[curcell,i])/nperms /sum(abs(rhat[curcell,])>cutp.g[curcell,i])
#			ncall.g[curcell,i]=sum(abs(rhat[curcell,])>cutp.g[curcell,i])
#		}
#		if(alternative == 'greater') {
#			fdr.g[curcell,i]=sum(rhatperm[,curcell,]>cutp.g[curcell,i])/nperms /sum(rhat[curcell,]>cutp.g[curcell,i])
#			ncall.g[curcell,i]=sum(rhat[curcell,]>cutp.g[curcell,i])
#		}
#		if(alternative == 'less') {
#			# [RG] BUG FIX: should be < - cutp.g[curcell,i]
##			fdr.g[curcell,i]=sum(rhatperm[,curcell,]< -cutp.g[curcell,i])/nperms /sum(rhat[curcell,]>cutp.g[curcell,i])
#			fdr.g[curcell,i]=sum(rhatperm[,curcell,]< -cutp.g[curcell,i])/nperms /sum(rhat[curcell,] < -cutp.g[curcell,i])
#			ncall.g[curcell,i]=sum(rhat[curcell,]< -cutp.g[curcell,i])
#		}		
#    }
#  }
#  
#  fdr.g =pmin(fdr.g,1)
#  
#  for (j in 1:numcell)	 {
#    fdr.g[j,]=make.monotone(fdr.g[j,])
#  }	

  FDR <- .csSAM_fdr(rhat, rhatperm, alternative)

  return (list(fdr.g=FDR$fdr.g, rhatperm = rhatperm, cutp.g = FDR$cutp.g, rhat = rhat,ncall.g = FDR$ncall.g, alternative))  
}


#' Computes FDRs in csSAM Models
#' 
#'  
#' @inheritParams fdrCsSAM0
#' @param alternative Type of test to conduct.
#' If "all" (default) then all althernative hypothesis are tested.
#' @param ... other arguments passed to [csSAMfit].
#' @param verbose Toggles verbose messages.
#' 
#' @export
fdrCsSAM <- function (G, cc, y, rhat, nperms, alternative = c('all', 'two.sided', 'greater', 'less'), ..., verbose = TRUE) {
  numgene <- ncol(G)
  numcell <- nrow(rhat)
  rhatperm <- array(dim = c(nperms,numcell,numgene))
  progress <- iterCount(n=nperms, title='', verbose = verbose == 1L)
  lapply(1:nperms, function(i, ...){
        progress()
        
        # permutation honouring blocks if any
        #ystar = y[sample_design(length(y), block = block)]
        ystar <- sample(y)
        rhatperm[i,,] <<- .csSAMfit(G = G, cc = cc, y = ystar, ..., verbose = max(verbose - 1L, 0L), stat.only = TRUE)
        NULL
      }, ...)
  progress(nperms, appendLF = FALSE)
  
  # compute fdr for all alternatives
  alternative <- match.arg(alternative)
  alternative <- unique(alternative)
  if( alternative == 'all' ) alternative <- c('two.sided', 'greater', 'less')
  FDR <- sapply(alternative, function(alt, ...){
        vmessage("Alternative '", alt,"' ... ", appendLF=FALSE)
        fdr <- .csSAM_fdr(..., alternative = alt)
        fdr$pvalue <- csSAM_pvalue(..., alternative = alt)
        vmessage('OK')
        fdr
      }
      , rhat = rhat, rhatperm = rhatperm
      , simplify = FALSE)
  
  return (list(FDR = FDR, rhatperm = rhatperm))
}


#' @useDynLib csSAM
.csSAM_fdr <- function(rhat, rhatperm, alternative, numcell = nrow(rhat), nperms = dim(rhatperm)[1L]){
    
    cutp.g=matrix(NA,nrow=numcell,ncol=100)
    
    for(j in 1:numcell)
        cutp.g[j,]=seq(0,max(abs(rhat[j,]), na.rm = TRUE),length=100)
    
    nperm_total <- dim(rhatperm)[1L]
    res <- .Call('csSAM_fdr', rhat, rhatperm, alternative, numcell, nperms, nperm_total, cutp.g, PACKAGE = 'csSAM')
    
    fdr.g <- res$fdr
    ncall.g <- res$ncall
    fdr.g =pmin(fdr.g, 1, na.rm = TRUE)
    
    for (j in 1:numcell)	 {
        fdr.g[j,]=make.monotone(fdr.g[j,])
    }	
    
    list(fdr.g=fdr.g, cutp.g = cutp.g, ncall.g = ncall.g, alternative = alternative)
}

# @param x statistics computed on original data (coef x feature)
# @param x statistics computed on permuted data (perms x coef x feature)
# @param alternative alternative hypothesis to test
permutation_pvalue <- function(x, xstar, alternative = c('two.sided', 'greater', 'less')){
  
  alternative <- match.arg(alternative)
  ncoef <- nrow(x)
  nfeatures <- ncol(x)
  nperms <- dim(xstar)[[1L]]
  pval <- matrix(NA_real_, ncoef, nfeatures, dimnames = list(rownames(x), colnames(x)))
  
  for (curcell in 1:ncoef) {
    if(alternative == 'two.sided') {
      np <- sweep(abs(xstar[,curcell,]), 2L, abs(x[curcell, ]), '-') >= 0
    }
    if(alternative == 'greater') {
      np <- sweep(xstar[,curcell,], 2L, x[curcell, ], '-') >= 0
    }
    if(alternative == 'less') {
      np <- sweep(xstar[,curcell,], 2L, x[curcell, ], '-') <= 0
    }	
    pval[curcell, ] <- (colSums(np) + 1) / (nperms + 1)
  }
  
  pval	
}

csSAM_pvalue <- function(rhat, rhatperm, ...) permutation_pvalue(rhat, rhatperm, ...)


