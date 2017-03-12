
#' @import utils stats Rcpp pkgmaker
NULL

#' Internal csSAM functions
#' 
#' These functions are not to be called by the user.
#' 
#' @author Shai Shen-Orr, Rob Tibshirani, Narasimhan Balasubramanian, David Wang
#' @keywords internal
#' @name csSAM-internal
NULL

#' csSAM package - Cell Type-Specific Statistical Analysis of Microarray
#' 
#' @name csSAM-package
#' @docType package
#' @details
#' Tissues are often made up of multiple cell-types. Each with its own
#' functional attributes and molecular signature.  Yet, the proportions of any
#' given cell-type in a sample can vary markedly. This results in a significant
#' loss of sensitivity in gene expression studies and great difficulty in
#' identifying the cellular source of any perturbations. Here we present a
#' statistical methodology (cell-type specific Significance Analysis of
#' Microarrays or csSAM) which, given microarray data from two groups of
#' biological samples and the relative cell-type frequencies of each sample,
#' estimates in a virtual manner the gene expression data for each cell-type at
#' a group level, and uses these to identify differentially expressed genes at a
#' cell-type specific level between groups.
#' 
#' The lower limit for the number of samples needed for deconvolving the
#' cell-specific expression of N cell-types is N+1. For a single colour array -
#' the result could be interpreted as the average expression level of a given gene
#' in a cell-type of that group. Multiplied by the frequency of a given cell-type
#' in an individual in the group, it is the amount contributed by that cell type
#' to the overall measured expression on the array.  
#' 
#' @section Key functions for this package:
#' \itemize{
#' \item \code{\link{csSamWrapper}} - Single wrapper function performs all
#' functionality. csfit: For deconvolving the average cell-type specific
#' expression for each cell-type in a given group.
#' \item \code{\link{csSAMfit}} - For fitting csSAM model on [Biobase::ExpressionSet] objects
#' with the convenience of a formula interface.
#' \item \code{\link{csSAM}} - For calculating the
#' contrast between every pair of cells being compared between the two
#' groups.
#' \item \code{\link{fdrCsSAM}} - Estimate the false discovery rate for each cell-type
#' specific comparison.
#' \item \code{\link{findSigGene}} - Identifies the list of differentially
#' expressed genes in a given cell-type at a given FDR cutoff.
#' \item \code{\link{plotCsSAM}} - Plots a fdr plot of their results.
#' }
#' 
#' Additional functions exists (e.g., \code{\link{runSAM}} and \code{\link{fdrSAM}} 
#' to contrast csSAM with the tissue heterogeneity ignorant SAM).
#' 
#' @bibliography ~/Documents/articles/library.bib
#' @cite Shen-Orr2010
#' @examples
#' 
#' library("csSAM")
#' ##
#' ## Generate random dataset
#' ##
#' set.seed(143)
#' k <- 5 # number of cell types 
#' ng <- 500 # number of genes
#' p <- 20 # number of samples
#' ndiff <- 100 # number of genes differentially expressed
#' 
#' # true cell-specific signatures
#' H1 <- matrix(rnorm(5*ng), ncol=ng)
#' H2 <- H1
#' # create differential expression for 3rd cell type
#' H2[3,1:ndiff] <- H2[3,1:ndiff] + 5 
#' 
#' # cell frequency matrix per sample
#' cc <- matrix(runif(p*k), ncol=k) 
#' cc <- t(scale(t(cc), center=FALSE, scale=rowSums(cc))) 
#' colnames(cc) <- paste('cellType', 1:ncol(cc), sep="")
#' 
#' # global expression matrix
#' G <- rbind(cc[1:10, ] %*% H1, cc[11:p, ] %*%H2 ) + matrix(rnorm(p*ng), ncol=ng)
#' # sample classes (2 groups) 
#' y <- gl(2, p/2)  
#' 
#' fileName = "Example File.pdf";
#' \dontshow{ on.exit(unlink(filename)) }
#' 
#' # Now run, either using the wrapper 
#' # NB: more permutations would be needed for real data
#' deconvResults = csSamWrapper(G, cc, y, nperms = 50, alternative = "two.sided"
#' 								, standardize = TRUE
#' 								, medianCenter = TRUE
#' 								, fileName = fileName)
#' 
#' # Or by calling each function independently: 
#' # this is useful if you want to perform only cell-specific expression 
#' # without differential expression.
#' 
#' # run csSAM
#' fit <- csSAMfit(t(G) ~ y | cc, nperms = 200, verbose = TRUE, .rng = 123)
#' top <- csTopTable(fit, alternvative = 'two.sided')
#' head(top[[1]])
#' 
#' # run SAM
#' tt.sam <- runSAM(G, y)
#' falseDiscovRSAM <- fdrSAM(G, y, tt.sam, nperms=200, alternative = 'two.sided')
#' 
#' plotCsSAM(fit, falseDiscovRSAM, alternative='two.sided', fileName)
#' 
NULL



