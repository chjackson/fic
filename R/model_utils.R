##' Form indicator matrix describing all possible submodels of a wide model
##'
##' @param wide A fitted model of standard R format, such that `terms(wide)` returns information about the terms of the model formula.  TODO do msm and flexsurv models support this?  msm probably not
##'
##' @param inds0 Narrow model indicators, in format described in \code{\link{fic}}
##'
##' @return A matrix in the format required by the `inds` argument of fic(), representing all possible submodels of the wide model.  The number of rows is the number of models, and the number of columns is the number of terms in the wide model.   The r,s entry of the matrix is a 1 if the rth submodel includes variable s, and 0 otherwise.   A factor is considered to be one term, i.e. there are not separate terms for each factor contrast.
##'
##' 
##' @export
all_inds <- function(wide, inds0=NULL){
    tt <- terms(wide)
    faclevs <- .getXlevels(tt, model.frame(wide))  # list of levels of all factors
    varnames <- names(attr(tt, "dataClasses"))[-1]
    contnames <- setdiff(varnames, names(faclevs)) # names of continuous variables
    contlevs <- Map(function(x)c(0,1), contnames)
    levs <- c(faclevs, contlevs)[varnames]
    dat <- do.call(expand.grid, levs)
    form <- delete.response(terms(wide))
    combs <- model.matrix(form, data=dat)[,-1,drop=FALSE]
    ## Keep only models which include the narrow model
    ininds0 <- if (!is.null(inds0)) apply(combs, 1, function(x) !any(x==0 & inds0==1)) else rep(TRUE, nrow(combs))
    combs[ininds0,]
}



##' Convert `newdata` to a design matrix excluding the intercept
##'
##' Numerics are allowed for character factor levels
##' 
getX <- function(newdata, wide){
    tt <- terms(wide)
    faclevs <- .getXlevels(tt, model.frame(wide)) # list of levels of all factors
    fn <- names(faclevs)
    newdata[fn] <- lapply(newdata[fn], as.character) ## convert numerics supplied for factor levels to characters
    X <- model.matrix(delete.response(tt), data=newdata, xlev=faclevs)[,-1,drop=FALSE]
    X
}
