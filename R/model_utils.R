##' Form indicator matrix describing all possible submodels of a wide model
##'
##' @param wide A fitted model of standard R format, such that `terms(wide)` returns information about the terms of the model formula.  Models outside standard R packages may not support this.
##'
##' @param inds0 Narrow model indicators, in format described in \code{\link{fic}}
##'
##' @return A matrix in the format required by the `inds` argument of fic(), representing all possible submodels of the wide model.  The number of rows is the number of models, and the number of columns is the number of parameters in the wide model.   The r,s entry of the matrix is a 1 if the rth submodel includes parameter s, and 0 otherwise.  If a factor is included (excluded) from the submodel, then all corresponding parameters are included (excluded) 
##'
##' 
##' @export
all_inds <- function(wide, inds0=NULL, intercept=TRUE){
    tt <- terms(wide)
    labs <- attr(tt, "term.labels")
    if (intercept) labs <- c("(Intercept)",labs)
    inds_short <- do.call(expand.grid, Map(function(x)c(0,1), labs))
    ass <- attr(model.matrix(wide), "assign")
    ass <- match(ass, unique(ass))
    combs <- inds_short[,ass,drop=FALSE]
    ininds0 <- if (!is.null(inds0)) apply(combs, 1, function(x) !any(x==0 & inds0==1)) else rep(TRUE, nrow(combs))
    combs[ininds0,]
}



##' Convert `newdata` to a design matrix excluding the intercept
##'
##' Numerics are allowed for character factor levels
##' 
newdata_to_X <- function(newdata, wide){
    tt <- terms(wide)
    faclevs <- .getXlevels(tt, model.frame(wide)) # list of levels of all factors
    fn <- names(faclevs)
    newdata[fn] <- lapply(newdata[fn], as.character) ## convert numerics supplied for factor levels to characters
    X <- model.matrix(delete.response(tt), data=newdata, xlev=faclevs)[,-1,drop=FALSE]
    X
}
