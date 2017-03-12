# Formula handling
# 
# Author: Renaud Gaujoux
###############################################################################

# internal wrapper
pData <- function(object){
  if( requireNamespace('Biobase') && isExpressionSet(object) ) Biobase::pData(object)
  else NULL
  
}

`%||%` <- function(a, b) if( is.null(a) ) b else a

#' @importFrom methods is
isExpressionSet <- function(x) is(x, 'ExpressionSet')
isExpressionMix <- function(x) is(x, 'ExpressionMix')

# @param object formula
# @param data cell type proportions matrix (cell types x samples)TODO: change this to a more standard variables in columns/data.frame input handling
#' @import Formula
csFormula <- function(object, data = NULL){
  
  x <- data
  # create multiple lhs/rhs formula object
  f <- f0 <- as.Formula(object)
  
  # check for response
  if( is.null(attr(f, 'lhs')) ){
    stop("Formula '", f,"' is invalid (no response term): target data must be specified as the left-hand/response term.")
  }
    
  # BASIS PART
  f.basis <- formula(f, rhs = 1)
  t <- terms(f.basis)
#   str(t)
  fenv <- environment(f.basis)
  # get target object
  vars <- attr(t, 'variables')
  ir <- attr(t, 'response')
  object <- eval(vars[[ir + 1L]], envir = fenv)
  
  # get factor as last covariate
  data <- NULL
  if( length(vars) > 2L ){
    tvar <- tail(vars, 1L)[[1L]]
    data <- if( isExpressionSet(object) ) eval(tvar, enclos = fenv, envir = pData(object))
        else eval(tvar, envir = fenv)
  }
  #
  
  # covariates
  covariates <- NULL
  coformula <- NULL
  if( length(vars) > 3L ){
    coformula <- update(formula(f, lhs = 0, rhs = 1), as.formula(sprintf("~ . - %s", as.character(tvar))))
    # use phenotypic data as formula data
    p.data <- if( isExpressionSet(object) ) pData(object)
    # add cell type proportions in formula data
    if( isExpressionMix(object) ) p.data <- cbind(p.data, as.data.frame(t(coef(object))))
    covariates <- model.matrix(coformula, data = p.data, na.action = NULL)
#        mf <- model.frame(coformula, data = p.data, na.action = NULL)
#        if( attr(mf, "dataClasses") ) covariates <- covariates[, -1L, drop = FALSE]
  }
  
  # load coefficient data (cell type proportions) and bind it with pheno data for
  # formula evaluation
  dfbind <- function(x, y){
    as.data.frame(append(x, y), row.names = rownames(x) %||% rownames(y))
  }
  coef_data <- try(coef(object), silent = TRUE)
  if( is(coef_data, 'try-error') ) coef_data <- NULL
  else if( !is.null(coef_data) ) coef_data <- as.data.frame(t(coef_data))
  f.data <- dfbind(coef_data, pData(object))
  
  # solve '.' into all proportions stored in LHS object as coefficients extracted with coef(lhs)
  rhs <- attr(f, 'rhs')
  if( rhs_is_dot <- (length(rhs) > 1L && identical(rhs[[2L]], quote(.))) ){
    varnames <- if( !is.null(coef_data) ) rownames(coef_data)
                else if( isExpressionSet(object) ){
                  # use all other variables
                  setdiff(colnames(pData(object)), attr(t, 'term.labels'))
                  
                }
    # expand formula 2nd RHS
    if( !is.null(varnames) ){
      attr(f, 'rhs')[[2L]] <- eval(parse(text = sprintf("quote(%s)", paste0(varnames, collapse = " + "))))
    }
  }
  
  # COEF PART
  # try to get x from the formula, otherwise ensure x is not missing (make it NULL)
  if( is.null(x) && length(rhs) > 1L ){
    xf <- model.frame(formula(f, lhs = 0, rhs = 2), data = f.data, na.action = NULL)
    x <- as.matrix(xf)
    # use actual cell type names for variable
    rn <- rownames(xf[[1L]]) %||% colnames(object)
    cn <- unlist(lapply(xf, colnames), use.names = FALSE)
    if( !length(cn) ) cn <- colnames(x)
    dimnames(x) <- list(rn, cn)
    
  } else if( !is.null(x) ) x <- t(x)
  ##

#    print(object)
#    str(data)
#    str(x)
#    stop()    
  list(formula = formula(f0), model = f, target = object, x = x, data = data, nvars = length(vars) - 2L, covariates = covariates, coformula = coformula)
}


