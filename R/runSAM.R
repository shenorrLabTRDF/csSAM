#' runSAM
#' 
#' A lightweight version of the SAM algorithm, only performs two group
#' comparison with equal deltas on each tail
#' 
#' 
#' @param G Matrix of gene expression, rows ordered in the same order at the
#' cell-frequency matrix (n by p, n samples, p genes)
#' @param y Numeric group association of each sample. Either 1 or 2.
#' @param s0.sam Input or computed value of SAM exchangeability factor. Default
#' is determined automatically
#' @param stand.r Median center and standardize arrays. Default is TRUE.
#' @param stat.only logical that indicates if only the SAM statistic should be 
#' returned, or a list containing a statistic matrix (SAM statistic, effect size and standard 
#' deviation) and the SAM exchangeability factor used in the computation.
#' @author Shai Shen-Orr, Rob Tibshirani, Narasimhan Balasubramanian, David Wang
#' @cite Shen-Orr2010
#' @export 
runSAM <-
function(G,y,s0.sam=NULL,stand.r=TRUE, stat.only = TRUE) {
  # BUGFIX: ttest.func works with genes in rows x samples in columns
  G <- t(G)
  if(stand.r == TRUE) {
	G=scale(G,center=apply(G, 2L, median),scale=F)
	}
  if(is.null(s0.sam)){
	s0.sam=quantile(ttest.func(G,y)$stat[, 'sd'],.5,na.rm=TRUE)
	}
  tt.sam=ttest.func(G,y,s0=s0.sam)
  if( stat.only ) return (tt.sam$stat[, 1L])
  else tt.sam
}

