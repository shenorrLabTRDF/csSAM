# Interface to Fit csSAM Models
# 
# Author: Renaud Gaujoux
###############################################################################

# TODO: move this function to pkgmaker or other generic package
.aggregate_list <- function(x){
  
  .tips <- function(x){
    unlist(sapply(names(x), function(n){
              if( is.list(x[[n]]) ) file.path(n, .tips(x[[n]])) else n
            }), use.names = FALSE)
  }
  
  .traverse <- function(x, name){
    res <- x
    lapply(name, function(n) res <<- res[[n]])
    res
  }
  
  t <- .tips(x[[1L]])
  tl <- strsplit(t, "/")
  res_flat <- lapply(tl, function(n){
        xl <- lapply(x, .traverse, n)
        fun <- if( !is.null(dim(xl[[1]])) ){
              if( nrow(xl[[1]]) == 1L ) rbind else cbind
            }else `c`
        do.call(fun, xl)
      })
  res_flat <- setNames(res_flat, t)
#      res_flat <- sapply(res_flat, function(x){
#                  if( is.null(dim(x)) && all(x == x[[1L]]) ) x[1L]
#                  else x
#              }, simplify = FALSE)
  
  # reshape into a nested list
  res <- list()
#    .nest <- function(n, np = NULL){
#        tip <- paste0(np, n, collapse = "/")
#        if( !length(n) ) return(res_flat[[tip]])
#        if( length(n) > 1 ){
#            res[[tip]] <<- .nest(n[-1], c(np, n[1]))
#        }else res[[n]] <<- res_flat[[tip]]
#    }
  
  .nest <- function(flat, children){
    if( !length(flat) ) return() 
    s <- split(seq_along(flat), sapply(children, head, 1L))
    res <- sapply(s, function(i){
          ch_i <- children[i]
          l1 <- lengths(ch_i) == 1L
          elmt <- flat[i[l1]]
          elmt <- c(elmt, .nest(flat[i[!l1]], lapply(ch_i[!l1], tail, -1L, simplify = FALSE)))
          if( length(elmt) == 1L ) elmt <- elmt[[1]]
          elmt
#       if( length(j) == 1L ) res_flat[[j]]
#       else .nest(j, lapply(children[j], tail, -1L, simplify = FALSE))
        }, simplify = FALSE)
    res
  }
  .nest(res_flat, tl)[names(x[[1]])]
}

# internal definition of basis() (original definition is in NMF)
basis <- function(x) x$basis
`basis<-` <- function(x, value){
  x$basis <- value
  x
}
basisfit <- function(x) x$basisfit

#' Fits a csSAM Model
#' 
#' 
#' @param x target matrix assayed on mixed samples.
#' @param ... arguments passed down to the relevant S3 method.
#' 
#' @return a `csSAMfit` object, that is a list with the following elements:
#'   *  
#'   * 
#'   *
#' @export
csSAMfit <- function(x, ...){
  UseMethod('csSAMfit')
}

#' @describeIn csSAMfit convenient formula interface that can take advantage of expression and 
#' sample annotation data stored in [Biobase::ExpressionSet] objects.
#' 
#' @export
csSAMfit.formula <- function(x, data = NULL, ...){
  # parse formula
  model <- csFormula(x, data = data)
  # call workhorse function
  csSAMfit.default(model$target, cc = t(model$x), y = model$data, covariates = model$covariates, ...)
  
}

#' @export
csSAMfit.Formula <- csSAMfit.formula

#' @describeIn csSAMfit workhorse method that takes model elements in separate
#' arguments.
#' 
#' @inheritParams csSamWrapper
#' @param covariates matrix of covariates to include in the model, e.g., for correcting 
#' for gender, age or batch effects.
#' @param fit fitting method, that dictates the type of model that is fitted:
#'   * auto: automatic selection based on input data and model specification (e.g., presence of covariates)
#'   * block: two linear models are fitted within each group of samples. Only works when the variable
#' of interest is a factor (categorical variable).
#'   * lm: full interaction linear model.
#'   * monovariate: model is fitted including one cell type at a time.
#' This could be a work around for when there are many more cell types than samples in a given group.
#' 
#' @param keep.all logical that indicates if the statistics of all the permutations should be also
#' returned in slot `rhatperm`.
#' @param verbose single logical or integer that controls verbosity levels.
#' @param .rng specification of the initial RNG settings, that are set using [rngtools::setRNG].
#' 
#' @import pkgmaker
#' @importFrom rngtools setRNG
#' @export
csSAMfit.default <- function(x, cc, y = NULL, covariates = NULL
    , nperms = 200, alternative = c('all', 'two.sided', 'greater', 'less')
    , fit = c('auto', 'block', 'lm', 'monovariate')
    , standardize=TRUE, medianCenter=TRUE
    , logRm =FALSE, logBase = 2, nonNeg=NULL
    , keep.all = FALSE, verbose = FALSE
    , .rng = NULL, ...){ 
  
  # check for extra arguments
  if( length(extra <- list(...)) )
    stop("Unused arguments: ", str_out(names(extra), Inf))
  
  # seed computation if requested
  if( !is.null(.rng) ){
    orng <- setRNG(.rng)
    on.exit( setRNG(orng) )
  }
  
  # set local verbose level
  if( !missing(verbose) || verbose ){
    ol <- lverbose(verbose)
    on.exit(lverbose(ol), add=TRUE)
  }
  
  
  ## map arguments to those of csSamWrapper
  if( is.null(cc) ) stop("Cannot fit csSAM model: missing proportion data.")
  
  Y <- x; x <- cc
  
  if( isExpressionSet(Y) ){
    if( !requireNamespace('Biobase') )
      stop("Missing dependency: package 'Biobase' is required to fit csSAM model on ExpressionSet objects.")
    Y <- Biobase::exprs(Y)
  }
  if( isExpressionMix(x) ) x <- coef(x)
  # data contains the definition of the effect of interest
  if( is.null(y) ) y <- factor(rep(1, ncol(Y)))
  if( !is.numeric(y) ){
    if( !is.factor(y) ) y <- factor(y, levels=unique(y))
    # remove absent levels
    y <- droplevels(y)
    
  } else stop("Invalid effect variable y: model does not handle variables of class [", class(y), ']')
  
  # transpose/copy for csSAM
  G <- t(Y); cc <- t(x)
  ##
  
  # show details of effect variable
  if( is.factor(y) ){
    # group composition
    n <- summary(addNA(y), maxsum=Inf)
    vmessage("Groups: ", str_out(n, Inf, use.names=TRUE, sep=' | '))
    
  } else{
    # range
    n <- summary(y)
    vmessage("Variable: ", str_out(n, Inf, use.names=TRUE, sep=' | '))
    
  }
  # show cell types
  vmessage("Cell type(s): ", if( is.null(colnames(cc)) ) ncol(cc) else str_out(colnames(cc), total = TRUE))
  
  # remove NA-labeled samples
  if( length(i_rm <- which(is.na(y))) ){
    wnote('Dropping samples with NA group/value: ', str_out(rownames(G)[i_rm]), ' [', length(i_rm), ']')
    y <- y[-i_rm]
    cc <- cc[-i_rm, , drop=FALSE]
    G <- G[-i_rm, , drop=FALSE]
    
  }
  
  res_object <- structure(list(basis = numeric(), coefficients = t(cc)), class = 'csSAMfit')
  
  # perform monovariate fit if requested
  fit <- match.arg(fit)
  vmessage("Fitting mode: ", fit)
  if( fit == 'monovariate' && ncol(cc) > 1L ){
    
    res <- matrix(NA, nrow(Y), ncol(cc), dimnames = list(rownames(Y), colnames(cc)))
    ca <- match.call()
    ca[[1L]] <- as.name('csSAMfit')
    ca[['x']] <- t(G)
    ca[['y']] <- y
    ca[['verbose']] <- max(verbose-1, 0)
    progress <- iterCount(n = ncol(cc), title='  Fitting monovariate models', extra = colnames(cc), verbose = verbose)
    bfit <- sapply(colnames(cc), function(ct){
          progress()
          ca[['cc']] <- t(cc[, ct, drop = FALSE])
          mono_fit <- eval(ca)
          res[, ct] <<- basis(mono_fit)
          mono_fit
          
        }, simplify = FALSE)
    progress(ncol(cc))
    
#    d <- dimnames(x)
    res_object$basisfit <- .aggregate_list(sapply(bfit, basisfit, simplify = FALSE))
    basis(res_object) <- res
#    dimnames(x) <- d
    
    return(res_object)
  }
  
  # remove samples that have missing data in any one of the cell types
  if( length(i_rm <- which(apply(cc, 1L, anyNA))) ){
    wnote('Dropping samples with missing proportions: ', str_out(rownames(G)[i_rm], total = TRUE))
    y <- y[-i_rm]
    cc <- cc[-i_rm, , drop=FALSE]
    G <- G[-i_rm, , drop=FALSE]
  }
  
  # show details after filtering
  if( is.factor(y) ){
    # group composition
    if( !identical(n, n2 <- summary(addNA(y), maxsum=Inf)) ) 
      vmessage("Groups (filtered): ", str_out(n2, Inf, use.names=TRUE, sep=' | '))
    
  }else{
    # range
    if( !identical(n, n2 <- summary(y)) )
      vmessage("Variable (filtered): ", str_out(n2, Inf, use.names=TRUE, sep=' | '))
    
  }
  n <- n2
  
  vmessage("Data (filtered): ", sprintf("%s features x %s samples", ncol(G), nrow(G)))
  
  # choose version based on inputs
  if( fit == 'auto' ){
    # covariates require version 2 (interaction model fit)
    fit <- 'block'
    if( !is.null(covariates) ){
      vmessage("Model has extra covariates: fitting lm interaction model")
      fit <- 'lm'
      
    } else if( is.numeric(y) ){
      vmessage("Model has numeric effect: fitting lm interaction model")
      fit <- 'lm'
      
    } else if( is.factor(y) && nlevels(y) > 2L ){
      vmessage("Model has factor effect with more than 2 levels: fitting lm interaction model")
      fit <- 'lm'
      
    }
    
  }
  version <- c(block = 1, lm = 2, monovariate = 1)[fit]
  
  # use factor encoding
#  if( is.factor(y) && version == 1 ) y <- as.integer(y)
  
  # detect negative values
  if( version == 1 && is.null(nonNeg) ){
    if( nonNeg <- (min(G, na.rm = TRUE) * min(cc, na.rm = TRUE) >= 0) ) vmessage("Input values: positives")
    else vmessage("Input values: mixed signs")
  }
  vmessage("Fitting model with ", if( isFALSE(nonNeg) ) "mixed sign" else "nonnegative", " effects")
  
#	numset = length(unique(y))
  numset <- nlevels(y)
  # check parameters
  if( numset > 2L ){
    if( version != 2 ){
      stop("csSAM - Cannot handle more than 2 groups of samples [", numset, "]")
    }else{
      vmessage("Model with more than 2 groups: switching to version 2")
      version <- 2
    }
    
  }
  if( nrow(G) != nrow(cc) ){
    stop("csSAM - Incompatible dimensions: number of samples in the"
        , " proportion matrix [", nrow(cc), "]"
        , " should match the number of samples in the target matrix [", nrow(G),"]")
  }
  #
  
  # fit csSAM model
  fitObject <- .csSAMfit(G, cc, y, covariates = covariates, version = version, stat.only = FALSE
      , standardize = standardize, medianCenter = medianCenter, nonNeg = nonNeg, logRm = logRm, logBase = logBase)
 
  # compute FDR
  if( nperms && is.null(fitObject$stat) ){
    warning("No cell type-specific differences were computed: FDR computation was skipped")
    nperms <- 0L
  }
  if( nperms ){
    alternative <- match.arg(alternative, several.ok = TRUE)
    vmessage("Computing FDR using ", nperms, " permutations ... ", appendLF=FALSE)
    elapsed <- system.time({
          fdr.csSAM <- fdrCsSAM(G,cc,y, fitObject$stat
              , nperms = nperms, alternative = alternative
              , standardize = standardize, medianCenter = medianCenter
              , logRm = logRm, logBase = logBase, nonNeg = nonNeg
              , verbose = max(verbose - 1L, 0), version = version, covariates = covariates
#              , block = block, rblock = rblock
          )
#          fdr.csSAM$alternative <- alternative
        })
    vmessage("OK")
    lmessage(2L, "Timing:")
    if( verbose >= 2L ) show(elapsed)
    
    # remove permuation results if not requested (to avoid memory issue)
    if( !keep.all ) fdr.csSAM$rhatperm <- NULL
  }
  
  
  # set result basis slot
  b <- t(fitObject$coefficients)
  rownames(b) <- rownames(Y)
  if( ncol(b) == ncol(cc) ){ # single interaction
    colnames(b) <- colnames(cc)
    basis(res_object) <- b
    
  }else{ # multi-interaction -> return an array of differences
    lev <- gsub(".*:([^:]+)$", "\\1", colnames(b))
    resA <- array(NA_real_, dim = c(nrow(b), ncol(cc), length(lev)), dimnames = c(rownames(b), colnames(cc), list(lev)))
    for(i in seq(lev)){
      resA[,,i] <- b[, grep(sprintf(":%s$", lev[i]), colnames(b))]
    }
    basis(res_object) <- resA
  }
  
  # attach fit result list
  res_object$basisfit <- c(fitObject, list(version = version, method = fit), if( nperms ) fdr.csSAM)
  res_object$call <- match.call()
  res_object$nperms <- nperms
  res_object$y <- y
  res_object$n <- ncol(G)
  
  # return
  res_object
}

# Internal workhorse function: assumes all parameters are in the correct format
# It is used when computing FDRs.
.csSAMfit <- function(G, cc, y, covariates = NULL, version = 1, standardize = TRUE, medianCenter = TRUE, nonNeg = TRUE, logRm = FALSE, logBase = 2
    , stat.only = TRUE, verbose = FALSE){
  
  # set local verbose level
  if( !missing(verbose) || verbose ){
    ol <- lverbose(verbose)
    on.exit(lverbose(ol), add=TRUE)
  }
  
  
  if( version == 2 ){
    vmessage("Fitting linear interaction model ... ", appendLF=FALSE)
    csSAMfit <- csSAM_single(cc, G, y, covariates = covariates, standardize = standardize, medianCenter = medianCenter, stat.only = stat.only)
    vmessage('OK')
    
  }else{
    vmessage("Fitting linear model within each group ... ", appendLF=FALSE)
#    sets <- split(1:nrow(G), y)
    n <- summary(y)
    deconv <- sapply(levels(y), function(l){
          i_set <- y %in% l
          cc_set <- cc[i_set, , drop=FALSE]
          G_set <- G[i_set, , drop=FALSE]
          if( nrow(G_set) <= ncol(cc_set) ){
            stop("csSAM - Insufficient sample size [", nrow(G_set), "]"
                , " in group '", l, "':"
                , " must be >= ", ncol(cc) + 1, " (i.e. number of cell types + 1).")
          }
          
          csfit(cc_set, G_set, logRm = logRm, logBase = logBase)
        }, simplify = FALSE)
    vmessage('OK')
    
    # early exit if only one group
    if( length(deconv) == 1L ){
      res_object <- list(stats = NULL)
      res_object$coefficients <- coef(deconv[[1L]]) 
      res_object$residuals <- residuals(deconv[[1L]])
      res_object$se <- deconv[[1L]]$se
      return(res_object)
    }
    
    #rhat <- array(dim = c(numcell,numgene))
    csSAMfit <- sapply(2:length(deconv), function(i){
          vmessage(sprintf("Computing csSAM model statistics [%s vs %s] ... ", names(n)[i], names(n)[1L]), appendLF=FALSE)
          csSAMfit <- csSAM(coef(deconv[[1L]]), deconv[[1L]]$se, n[1L]
              , coef(deconv[[i]]), deconv[[i]]$se, n[i]
              , standardize = standardize, medianCenter = medianCenter, nonNeg = nonNeg, stat.only = stat.only)
#          colnames(csSAMfit$stat) <- colnames(G)
          vmessage("OK")
          csSAMfit
        }, simplify = FALSE)
    if( nlevels(y) == 2L ) csSAMfit <- csSAMfit[[1L]]
    if( !stat.only ) csSAMfit$csfit <- deconv
  }
  
  # return result
  csSAMfit
}
