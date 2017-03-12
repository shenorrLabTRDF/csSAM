# Cell-specific statistics
# 
# Author: Renaud Gaujoux
# Created: 11 Dec 2012
###############################################################################


#' Compute Cell-Specific Statistics
#' 
#' \code{csTopTable} is a generic function that returns cell type-level statistics 
#' on biological features, e.g., cell type-specific expression.
#' 
#' @param x data object, typically obtained by fitting a model that computes 
#' cell type-specific statistics, e.g., [csSAMfit].
#' @param ... extra parameters passed to specific methods.
#' 
#' @export
csTopTable <- function(x, ...){
    UseMethod('csTopTable')
}

# returns margin index or names if defined
dimidx <- function(x, margin){
    d <- dimnames(x)[[margin]]
    if( is.null(d) ) d <- seq(1, dim(x)[[margin]]) 
    d
}

# returns the number of decimals in a real number
digits <- function(x){
  res <- rep(NA_integer_, length(x))
  res[grep('[.,]', x, invert = TRUE)] <- 0L
  res[is.na(res)] <- nchar(gsub(".*[.,]([^.,]+)$", "\\1", x[is.na(res)]))
  res
}


#' @param n maximum number of top features to retain -- in each cell type.
#' If `NULL` or `Inf` then all selected features are returned.
#' Otherwise the first `n` feature in each cell type are returned, **after**
#' ordering according to arguments `sort.by` and `decreasing`.
#' @param sort.by variable to use to order the features.
#' @param select specifies the statistic to include in the result.
#' @param subset indicates the subset of tissue/cell type to include in the result.
#' Default is to return only the cell types for which a `threshold` value was 
#' provided. If `NULL`, as default, then no subsetting is performed.
#' @param decreasing logical that controls the feature ordering based on `sort.by`.
#' @param threshold numeric vector that indicates the threshold applied to select 
#' features based on their values of `sort.by`.
#' Thresholds can be defined for specific cell types by using naming the elements 
#' of `threshold` with corresponding cell type names.
#' The default threshold that is applied to all cell types is taken from the first
#' unnamed element.
#' 
#' Depending on `decreasing`, thresholds are interpreted differently:
#' 
#'   * `decreasing = FALSE` keeps features with `sort.by` "<= x"
#'   * `decreasing = TRUE` keeps features with `sort.by` ">= x"
#' 
#' @param annotation Bioconductor annotation package to use to annotate 
#' the features (rownames) using [xbioc::geneInfo].
#' @param merge logical that indicates if the result should be a 
#' a list of `data.frame` objects (`merge = FALSE`), or a single 
#' `data.frame` object with columns "Feature" and "Cell type" containing
#' the feature identifiers and cell type respectively.
#' @param nodups a logical that indicates if features appearing in multiple 
#' cell types should be removed.
#' This is useful in combination with `threshold` to extract features that are 
#' found significant in a single cell type.
#' 
#' @return a list of `data.frame` objects if `merge = FALSE`, or a single 
#' `data.frame` object if `merge = TRUE`.
#' 
#' If `merge=TRUE` and `select` is of length one, then the result is a `matrix` object
#' (feature x cell type) containing the selected value.
#' 
#' @importFrom plyr ldply
#' @export
#' @rdname csTopTable
csTopTable.default <- function(x, n=100L, sort.by = 1L, decreasing=FALSE
                                , select = NULL, subset = names(threshold)
                                , threshold = NULL, annotation = NULL
                                , ...
                                , merge = FALSE, nodups = FALSE){
        
    if( !length(x) ) return(x)
    
    # quick disable any filtering/ordering
    if( is_NA(n) ){
        decreasing <- NA
        n <- NULL
    }
    # quick merge and return matrix of statistics
    if( is.null(select) && is.character(merge) ){
        select <- merge
        merge <- TRUE
    }
    
    # comparison operator
    th_cmp <- match.fun(ifelse(isTRUE(decreasing), '>=', '<='))
    # add names if necessary
    if( is.null(names(x)) ) names(x) <- as.character(seq_along(x))
    # expand threshold if necessary
    if( length(threshold) == 1L && is.null(names(threshold)) ){
      threshold <- setNames(rep(threshold, length(x)), names(x))
    }
    # define default threshold from first unlabeled threshold
    if( length(default_threshold_i <- which(names(threshold) %in% '')) ){
      threshold[setdiff(names(x), names(threshold))] <- threshold[default_threshold_i]
      threshold <- threshold[!names(threshold) %in% '']
    }
   
    # select specific cell types
    if( !is.null(subset) ){
      # handle ct:effect format
      if( all(grepl(':', names(x))) )
        i_subset <- na.omit(unique(pmatch(c(subset, paste0(subset, ':')), names(x))))
      x <- x[i_subset]
    }
    
    top <- sapply(names(x), function(ct){
            ct_res <- x[[ct]]
            if( is.matrix(ct_res) ) ct_res <- as.data.frame(ct_res)
            
            # split cell type into cell type and effect names
            ct_name <- ct
            if( grepl(":", ct) ){
              ct_res <- cbind(ct_res, Level = sub(".*:([^:]+)$", "\\1", ct), stringsAsFactors = FALSE)
              ct_name <- sub("(.*):[^:]+$", "\\1", ct)
              
            }
            ct_res <- cbind(ct_res, Cell.type = ct_name, stringsAsFactors = FALSE)
            
            # add rownames if necessary
            if( is.null(rownames(ct_res)) ) rownames(ct_res) <- 1:nrow(ct_res)
            
            # apply threshold
            th <- threshold[ct]
            if( is_NA(th) ) th <- threshold[ct_name]
            if( length(threshold) && !is_NA(th) ){
              # use which() to handle NA values
              ct_res <- ct_res[which(th_cmp(round(ct_res[, sort.by], digits(th)), th)), , drop = FALSE]
              # exit if nothing is left
              if( !nrow(ct_res) ) return(ct_res)
              
            }
            
            # order
            if( !is_NA(decreasing) ) ct_res <- ct_res[order(ct_res[, sort.by], decreasing = decreasing), , drop = FALSE]
            
            # truncate
            if( !is.null(n) && is.finite(n) ){
              ct_res <- head(ct_res, n)
              # exit if nothing is left
              if( !nrow(ct_res) ) return(ct_res)
              
            }
              
            # annotate
            if( !is.null(annotation) && requireNamespace('xbioc') ) 
              ct_res <- cbind(xbioc::geneInfo(rownames(ct_res), annotation = annotation, ...), ct_res)
            
            # return
            ct_res
          }, simplify = FALSE)
    
    
     # warns for unsused thresholds
    if( length(bad_names <- setdiff(names(threshold), c(names(x), unique(unlist(sapply(top, function(d) d[, 'Cell.type'])))))) ){
      warning("Unused thresholds for cell types not found in data: ", str_out(bad_names, total = TRUE))
    }
          
    
    # select specific column
    if( !is.null(select) ){
        simply <- length(select) == 1L
        top <- sapply(top, function(x) x[, select, drop = simply], simplify = simply)         
    }
        
    if( is.list(top) && merge ){
        all0 <- as.data.frame(top[[1]], stringsAsFactors = FALSE)
        all_data <- ldply(names(top)[sapply(top, nrow) > 0], function(ct){
                stats <- top[[ct]]
                res <- data.frame(Feature = rownames(stats), `Effect` = ct)
                cbind(res, stats)
              })
        top <- if( nrow(all_data) ) all_data else all0
        
        # only if there is anyting left
        if( nrow(top) ){
          # reorder
          if( !is_NA(decreasing) ){
              if( is.numeric(sort.by) ) sort.by <- sort.by + 2L
              top <- top[order(top[, sort.by], decreasing = decreasing), , drop = FALSE]
          }
          
          # only keep exclusive genes if requested  
          if( nodups ){
            dups <- top$Feature[duplicated(top$Feature)]
            top <- top[!top$Feature %in% dups, ]
          }
        }
  
    }
    
    # reorder columns
    .order <- function(v) order(match(v, c('Feature', 'Effect', 'Cell.type', 'Level')))
    if( !is.null(ncol(top)) ) top <- top[, .order(colnames(top)), drop = FALSE]
    else top <- sapply(top, function(t) t[, .order(colnames(t)), drop = FALSE], simplify = FALSE)
    
    # attach annotation package to result
    attr(top, 'annotation') <- annotation
    class(top) <- c('csTopTable', class(top))
    
    top
}

#' @export
csTopTable.matrix <- function(x, ...){
    
    # convert into a list
    cts <- dimidx(x, 2L)
    xl <- sapply(cts, function(i){ x[, i, drop = FALSE] }, simplify=FALSE)
    csTopTable(xl, ...)
}

#' @export
#' @rdname csTopTable
csTopTable.array <- function(x, ...){
    
    if( length(dim(x)) == 2L ) csTopTable(x[,], ...)
    else{
        cmp <- dimidx(x, 1L)
        sapply(cmp, function(i, ...){
                    csTopTable(x[i,,], ...)
                }, ..., simplify=FALSE)
    }
}

#' @export
#' @rdname csTopTable
csTopTable.character <- function(x, ...){
    cl <- match.call()
    cl[['x']] <- NULL
    cl[[1L]] <- as.name(paste0('csTopTable.', x))
    e <- parent.frame()
    eval(cl, e)
}

#' @describeIn csTopTable returns, for each feature, the false discovery 
#' rates of differential expression between groups of samples within each cell type, 
#' as computed by \code{\link[=csSAM]{fdrCsSAM}} when running csSAM.
#' These are returned as a list, whith one element per cell type.
#' 
#' @param alternative Alternative hypothesis to consider.
#' Alternative "min" returns, for each features, the alternative achieving the minimum 
#' FDR.
#' 
#' @seealso [fdrCsSAM]
#' 
#' @export
csTopTable.csSAMfit <- function(x, alternative = c('two.sided', 'greater', 'less', 'min'), sort.by = 'FDR', decreasing = FALSE, ...){
  
  # extract fitted model
  fit <- x$basisfit %||% x
  
  if( !is.list(fit) ) stop("Invalid input data: must be a list.")
  if( is.null(rhat <- fit$stat) ){
    stop("Cannot compute top table: csSAM fit does not contain group comparison data.\n  Was csSAM run with a group variable?")
  }
  
  # extract baseline levels
  if( !is.null(fit$fit) ){
    full_coef <- coef(fit$fit)
    base_coef <- setdiff(rownames(full_coef), rownames(rhat))
  }else{
    full_coef <- coef(fit$csfit[[1]])
    base_coef <- rownames(full_coef)
  }
  #
  
  res <- sapply(rownames(rhat), function(ct){
        
        # get corresponding baseline level
        i_base <- which( !is.na(pmatch(base_coef, ct)) )
        if( length(i_base) != 1L && !is.na(i_base) ) stop("Could not extract baseline coefficient for cell type '", ct, "'")
        i_base <- base_coef[i_base]
        #
        
        cbind(`|t|` = abs(rhat[ct, ]), t = rhat[ct, ]
            , Baseline = full_coef[i_base, ], Delta = coef(fit)[ct, ]
            , se = fit$se[ct, ])
      }, simplify = FALSE)
  
  if( !is.null(fit$FDR) ){
    # vmessage("Finding signature genes ... ", appendLF=FALSE)
    if( !missing(alternative) ) alternative <- match.arg(alternative, c('min', names(fit$FDR)))
    else alternative <- names(fit$FDR)[1L]
    
    .compute_fdr <- function(alternative){
      alt_code <- match(alternative, c('less', 'two.sided', 'greater')) - 2L
      csSAMData <- fit$FDR[[alternative]]
      fdr <- t(findSigGene(rhat = rhat, csSAMData = csSAMData, alternative = alternative))
      rownames(fdr) <- colnames(rhat)
      colnames(fdr) <- rownames(rhat)
      
      pvalue <- t(csSAMData$pvalue)
      sapply(colnames(fdr), function(ct){
            cbind(FDR = fdr[, ct], P.Value = pvalue[, ct], Alternative = alt_code, res[[ct]])
          }, simplify = FALSE)
    }
    
    if( alternative == 'min'){ # select min FDR across alternatives
      fdr <- sapply(names(fit$FDR), .compute_fdr, simplify = FALSE)
      res <- sapply(names(fdr[[1L]]), function(ct){
            res <- fdr[[1L]][[ct]]
            x <- sapply(fdr, function(x) x[[ct]][, 'FDR'])
            j <- max.col(-x, 'first')
            sapply(seq_along(j), function(i){
                  res[i, ] <<- fdr[[j[i]]][[ct]][i, ]
                })
            res
          }, simplify = FALSE)
    }else res <- .compute_fdr(alternative)
    
  }else if( missing(sort.by) ){ # sort by absolute effect t-statsitic
    sort.by <- '|t|'
    if( missing(decreasing) ) decreasing <- TRUE
    
  }
  # vmessage('OK')  
  
  top <- csTopTable(res, ..., sort.by = sort.by, decreasing = decreasing)
  class(top) <- c('csTopTable_csSAMfit', class(top))
  top
}
