#' csSAM
#' 
#' Computes the constrast between groups for the deconvolved cell-specific
#' expression for each cell-type
#' 
#' 
#' @param ghat1 Expression matrix of deconvolved cell-specific gene expression
#' estimates for group 1.
#' @param se1 Standard error group 1
#' @param n1 Group 1 size
#' @param ghat2 Expression matrix of deconvolved cell-specific gene expression
#' estimates for group 2.
#' @param se2 Standard error group 2
#' @param n2 Group 2 size
#' @param standardize Standardize contrast values, into SAM shrinked statistic.
#' @param medianCenter Median center rhat distributions for each cell-type
#' @param nonNeg Negative values not allowed such as in a single channel
#' microarray. Zero them if negative (a conervative option)
#' @param stat.only logical that indicates if the returned object should consists in the 
#' matrix of cell type-specific SAM statistic only, or contain more detailed data about 
#' the fit (see section \emph{Value}). 
#' 
#' 
#' @return A matrix object that contains the -- standardised -- contrast of the average
#' cell-specific expression profile of the two groups, per cell-type (Size k by
#' g where k is the number of cells and g is the number of genes).
#' 
#' If \code{stat.only=FALSE}, then result is a list with the following elements, where
#' matrix elements have cell types in rows and genes in columns:
#' \describe{
#' \item{coefficients}{a matrix of cell-specific contrasts between the two groups}
#' \item{se}{a matrix of cell-specific standard errors}
#' \item{s0}{a matrix of cell-specific standard error SAM shrinkage}
#' \item{stat}{a matrix of -- possibly median-centered -- cell-specific SAM shrinked t-test statistics}
#' \item{stat.med}{a vector containing the medians of each cell-specific SAM shrinked t-test statistics}
#' } 
#' 
#' @author Shai Shen-Orr, Rob Tibshirani, Narasimhan Balasubramanian, David Wang
#' @cite Shen-Orr2010
#' @export
csSAM <-
function(ghat1, se1, n1, ghat2, se2, n2, standardize,
                  medianCenter=TRUE, nonNeg=FALSE, stat.only = TRUE) {
  numcell=nrow(ghat1)
  if (nonNeg) {
    ghat1 <- pmax(ghat1, 0)	
    ghat2 <- pmax(ghat2, 0)
  }
  rhat = ghat2 - ghat1
  se=((n1*se1^2+n2*se2^2)/(n1+n2))^(1/2)
  csSAM_stats(rhat, se, standardize = standardize, medianCenter = medianCenter, stat.only = stat.only)
  
}

# Fit csSAM model as a single interaction model
csSAM_single <- function(cc, G, gr, covariates = NULL, standardize = TRUE, medianCenter = TRUE, stat.only = TRUE){
  
  # fit single regression model
  
  .design <- function(cc, gr, covariates = NULL){
    D <- model.matrix(~ cc + gr:cc)
    
    # rename coloumns to match cell types
    xsub <- ''
    if( !is.matrix(cc) ) deparse(substitute(cc))
    else if( ncol(cc) == 1L ) xsub <- colnames(cc)
    colnames(D) <- gsub("^cc([^:]*):gr", paste0(xsub, "\\1:"), colnames(D))
    colnames(D) <- gsub("^cc", xsub, colnames(D))
    
    # compute number of interaction coefficients: all - intercept - number of cell types 
    interact_i <- ncol(D) - 1 - ncol(cc)
    # handles extra covariates
    if( is.null(covariates) ) D <- D[, -1L, drop = FALSE]
    else{
#        str(covariates)
      D <- cbind(covariates, D[, -1L, drop = FALSE])
    }
    idiff <- seq(ncol(D) - interact_i + 1, ncol(D))
    list(D = D, idiff = idiff)
  }
  
  design <- .design(cc, gr, covariates)
  D <- design$D
  idiff <- design$idiff
  
#    print(colnames(D))
#    print(colnames(D)[idiff])
#    print(head(D, 20))
#    print(colnames(D))
  
  fit1 <- lsfit(D, G, intercept=FALSE) 
  se1 <- ls.diag(fit1)$std.err
  res <- structure(list(ghat = coef(fit1), residuals = residuals(fit1), se = se1, n = nrow(G))
      , class = 'csfit')
  # fix dimnames
  dimnames(res$residuals) <- dimnames(G)
  colnames(res$ghat) <- colnames(res$se) <- colnames(G)
  
  # compute SAM statistic for the coefficients of the cell type differences
  res <- csSAM_stats(res$ghat[idiff, , drop = FALSE], res$se[idiff, , drop = FALSE], standardize = standardize, medianCenter = medianCenter, stat.only = stat.only)
  
  if( !stat.only ){
    res$fit <- fit1
  }
  res
}

csSAM_stats <- function(rhat, se, standardize, medianCenter=TRUE, stat.only = TRUE) {
  
  numcell <- nrow(rhat)
  if( !stat.only ){
    res <- structure(list(coefficients = rhat, se = se, standardize = standardize), class = 'csSAM')
  }
  
  ##if the data is to be standardized
  if ( !isFALSE(standardize) ) {
    if( isTRUE(standardize) ) standardize <- 'shrink'
    
    if( !stat.only ){
      res$standardize <- standardize
    } 
    s0.r <- switch(standardize, 
        shrink = apply(se,1,quantile, .5,na.rm=TRUE)
        , 'std' = rep(0, numcell)
        , stop("Invalid standardize method: available methods are 'shrink' and 'std'."))
    
    for (i in 1:numcell) {
      rhat[i,]=rhat[i,]/(se[i,]+s0.r[i])
    }
    
    if( !stat.only ){# store std error and shrinkage
      res$s0 <- s0.r
    }
  }
  
  ##if rhat is to be median-centered
  if (medianCenter) {
    rmed <- apply(rhat, 1L, median, na.rm = TRUE)
    if( !stat.only ){
      res$stat.med <- rmed
    }
    rhat <- sweep(rhat, 1L, rmed, '-')
  }
  
  if( !stat.only ){
    res$stat <- rhat
  }
  
  if( !stat.only ) return(res)
  else return (rhat)
}

