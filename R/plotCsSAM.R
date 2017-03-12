#' plotcsSAM
#' 
#' Plots the # of genes called significnat at a given false disocvery rate for
#' the SAM (heterogenous tissue) comparison, and for each of the contrasted
#' cell-types using csSAM
#' 
#' 
#' @param csfit csSAM fitted model as returned by [csSAMfit].
#' @param sam List object output of the fdrSAM function
#' @param alternative Type of test conducted. Will appear in plot title.
#' @param fileName Name of output pdf file.
#' @param ... other arguments passed to [plot.csSAMfit].
#' @param height,width output height and width
#' @author Shai Shen-Orr, Rob Tibshirani, Narasimhan Balasubramanian, David Wang
#' @cite Shen-Orr2010
#' 
#' @importFrom grDevices dev.off pdf
#' @importFrom graphics plot par title
#' @export
plotCsSAM <- function(csfit, sam, alternative = 'two.sided', fileName, ..., height = 8, width = 10){
  
  
  fileName <- paste0(sub("\\.pdf$", "", fileName), '.pdf')
  on.exit( unlink(fileName) )
  pdf(file = fileName, height = height, width = width)

  # SAM
  n <- sam$ncall.sam
  fdr <- sam$fdr.sam[n>0]
  n <- n[n>0]
  # complete data if necessary
  if( min(n) != 1L ){
    n <- c(n, 1L)
    fdr <- c(fdr, min(fdr, na.rm = TRUE))
  }
  plot(n, fdr, xlab="## called", ylab="FDR", type="l",log="x",ylim=c(0,1))
  title(paste("SAM", alternative))
  # csSAM
  p <- plot(csfit, alternative = alternative, ...) + ggtitle('csSAM')
  print(p)
  
  on.exit()
  dev.off()
  # return path
  invisible(fileName)
  
}

#' Plotting csSAM FDR
#' 
#' @param x a csSAM fit as returned by [csSAMfit]
#' @param  subset character or index vector that specifies the subset of
#' cell types to include in the plot.
#' @param alternative test alternative hypothesis for which FDRs are shown.
#' Multiple alternatives can be plotted at once ("two.sided", "greater", "less").
#' Use `alternative = "all"` to plot all of them.
#' @param by specifies the facetting criteria.
#' Default is produces one panel per cell type, each showing FDR for the selected alternative(s).
#' @param xlab,ylab axis labels
#' @param ylim y-axis range
#' @param ... unused trailing arguments for S3 conformity.
#' An error will be thrown if any unsused argument are passed to the function.
#' 
#' @export
#' @import ggplot2
#' @importFrom scales comma
plot.csSAMfit <- function(x, subset = NULL, alternative = 'two.sided', by = c('cell-type', 'alternative')
                          , xlab='# called', ylab= 'FDR', ylim=c(0,1), ...){
  
  if( length(extra <- list(...)) )
    stop("Unused arguments: ", str_out(names(extra), Inf))
  
  # extract fitted model
  fit <- basisfit(x)
  
  if( !is.list(fit) ) stop("Invalid input data: must be a list.")
  if( is.null(rhat <- fit$stat) ){
    warning("Nothing to plot: data does not contain any group comparison data [was csSAMfit run with a group variable in argument `data`?]")
    return(invisible())
  }
  
  ncell <- nrow(rhat)
  cellnames <- rownames(rhat)
  if( is.null(cellnames) ){
    if( is.character(subset) ){
      if( length(subset) > ncell ){
        stop("Invalid number of cell types (", length(subset), "):"
            , " should have at most the same number of cell type data in `x` (", ncell, ")")
      }
      cellnames <- subset
    }else cellnames <- as.character(seq(ncell))
  }
  if( is.null(subset) ) subset <- seq(ncell)
  if( is.character(subset) ){ # partial match the types
    subset <- pmatch(subset, cellnames)
    subset <- subset[!is.na(subset)]
    if( !length(subset) )
      stop("None of the cell types ", str_out(subset), " matched [available: ", str_out(cellnames),"]")
  }
  

  if( is.null(FDR <- fit$FDR) )
    stop("Could not find FDR data in fit object: check if csSAM was run with group data and nperms > 0.")
#            stop("Could not find FDR data in fit object: provide reference data, or, alternatively, run csSAM with group data and nperms > 0.")
  
  alt_set <- names(FDR)
  if( !'all' %in% alternative ){
    ialt <- charmatch(alternative, alt_set, nomatch = 0L)
    FDR <- FDR[ialt]
  }
  if( !length(FDR) ){
    stop("No FDR data found for alternative ", str_out(alternative, Inf)
        , "\n  Available hypothesis: ", str_out(alt_set, Inf))   
  }
  
  df <- ldply(setNames(subset, cellnames[subset]), function(i){                
        df <- ldply(FDR, function(x){
              n <- ncol(x$ncall.g)
              j <- which(x$ncall.g[i, ]<=0)
              if( length(j) != n ){
                nc <- rev(x$ncall.g[i, -j])
                ifirst <- !duplicated(nc)
                val <- rev(x$fdr.g[i, -j])[ifirst]
                nc <- nc[ifirst]
                val <- as.numeric(rbind(val, c(val[-1], tail(val, 1L))))
                nc <- as.numeric(rbind(nc, nc))
                
                
                # complete plot with data for first called
                if( nc[1L] != 1 ){
                  nc <- c(1, nc)
                  val <- c(val[1], val) 
                }
                
                data.frame(x = nc, y = val)
              }else
                data.frame(x = seq(1, ncol(rhat), length.out = n), y = rev(x$fdr.g[i,]))
              
            }, .id = 'Alternative')
        df
      }, .id = 'Cell.type')
  
  # order cell type levels as input columns
  df[['Cell.type']] <- factor(df[['Cell.type']], levels = cellnames)
  
  
  # compute breaks in log scale
#  breaks <- .log10_break(df$x)
  
  by <- match.arg(by)
  if( by == 'alternative' ){
    by.colour <- 'Cell.type'
    by.facet <- 'Alternative'
  }else if( by == 'cell-type' ){
    by.colour <- 'Alternative'
    by.facet <- 'Cell.type'
  }
  
  aes_specs <- aes_string(colour = by.colour) 
  
  label_pretty <- as_labeller(function(labels){
    x <- gsub('.', ' ', labels, fixed = TRUE)
    if( any(w <- grepl(' ', x)) ) x[w] <- paste0(toupper(substr(x[w], 1, 1)), substring(x[w], 2))
    x
  })
  # create ggplot plot
  p <- ggplot(df, aes_string(x = 'x', y = 'y')) + 
      geom_line(aes_specs) +
      scale_y_continuous(ylab, limits=ylim) + 
#      scale_x_continuous(xlab, breaks = breaks) +
      scale_x_log10(labels = scales::comma) + 
      theme_bw() +
      theme(axis.text.x = element_text(angle = 50, hjust = 1)) +
      facet_grid(paste0('~ ', by.facet), labeller = label_pretty)  
  
  ltitle <- label_pretty(by.colour)
  p <- p + labs(colour = ltitle, linetype = ltitle)
  
  
  # return plot object
  p
}
