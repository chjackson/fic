##' Convert 0/1 indicators for inclusion of regression terms to
##' indicators for inclusion of parameters
##'
##' If a regression term is a factor, then the 0 or 1 indicating its
##' inclusion/exclusion is replicated to a length given by the number
##' of factor levels minus 1, that is, the number of parameters
##' pertaining to that factor.
##'
##' This function only works for classes of models for which the
##' \code{model.matrix} function is understood and returns objects
##' with an \code{"assign"} attribute. This includes all the
##' commonly-used models in base R.
##' 
##' @inheritParams fic
##' 
##' @export
expand_inds <- function(inds, wide){
    if (!(is.vector(inds) || is.matrix(inds) || is.data.frame(inds)))
        stop("`inds` must be a vector, matrix or data frame")
    if (is.data.frame(inds))
        inds <- as.matrix(inds)
    if (!is.numeric(inds)) stop("`inds` must contain only numeric elements")
    ass0 <- attr(model.matrix(wide), "assign")
    ass <- match(ass0, unique(ass0))
    nterms <- length(unique(ass))
    has_icpt <- (0 %in% ass0)
    icpt_msg <- if (has_icpt) ", including the intercept" else ""
    if(is.vector(inds)) {
        if (length(inds) != nterms){
            stop(sprintf("`inds` of length %s, but %s terms in model%s", length(inds), nterms, icpt_msg))
        }
        inds <- matrix(inds, nrow=1)
    } else {
        if (ncol(inds) != nterms){
            stop(sprintf("`inds` has %s columns, but %s terms in model%s", length(inds), nterms, icpt_msg))
        }
    }
    inds[,ass,drop=FALSE]
}



##' Form indicator matrix describing all possible submodels of a wide model
##'
##' @param wide A fitted model of standard R format, such that `terms(wide)` returns information about the terms of the model formula.  Models outside standard R packages may not support this.
##'
##' @param inds0 Narrow model indicators, in format described in \code{\link{fic}}
##'
##' @param intercept Is a regression intercept included in the indicators?  Should be \code{TRUE} for standard fully parametric regression models, and \code{FALSE} for Cox regression.
##'
##' @return A matrix in the format required by the `inds` argument of fic(), representing all possible submodels of the wide model.  The number of rows is the number of models, and the number of columns is the number of parameters in the wide model.   The r,s entry of the matrix is a 1 if the rth submodel includes parameter s, and 0 otherwise.  If a factor is included (excluded) from the submodel, then all corresponding parameters are included (excluded) 
##'
##' 
##' @export
all_inds <- function(wide, inds0=NULL, intercept=TRUE){
    tt <- terms(wide)
    labs <- attr(tt, "term.labels")
    if (inherits(wide, "coxph"))
        intercept <- FALSE
    if (intercept) labs <- c("(Intercept)",labs)
    inds_short <- do.call(expand.grid, Map(function(x)c(0,1), labs))
    ass <- attr(model.matrix(wide), "assign")
    ass <- match(ass, unique(ass))
    combs <- inds_short[,ass,drop=FALSE]
    rownames(combs) <- apply(combs, 1, paste, collapse="")
    ininds0 <- if (!is.null(inds0)) apply(combs, 1, function(x) !any(x==0 & inds0==1)) else rep(TRUE, nrow(combs))
    combs[ininds0,]
}



##' Convert data frame of covariate values to a design matrix, excluding the intercept
##'
##' @inheritParams
##' 
##' @param newdata Data frame where each row is a vector of covariate values defining an alternative focus quantity.
##'
##' @param wide Wide model which includes these covariates.
##'
##' @return "Design" matrix of covariate values defining alternative focuses, with factors expanded to their contrasts
##'
##' @details Numerics are allowed for character factor levels
##'
##' @export
newdata_to_X <- function(newdata, wide){
    tt <- terms(wide)
    faclevs <- .getXlevels(tt, model.frame(wide)) # list of levels of all factors
    fn <- names(faclevs)
    newdata[fn] <- lapply(newdata[fn], as.character) ## convert numerics supplied for factor levels to characters
    X <- model.matrix(delete.response(tt), data=newdata, xlev=faclevs)[,-1,drop=FALSE]
    X
}


##' Fit submodels of a wide model, defined by a matrix of indicators
##' for inclusion of covariates.
##' 
##' Fit the submodels of a wide model \code{wide} which are defined by
##' \code{inds}.  This can only be used for covariate selection
##' problems, where the submodels contain different subsets of
##' covariates.
##'
##' Requires \code{wide} to have a component named
##' \code{call} giving the function call used to produce \code{wide}.
##' This call should include a \code{formula} component, which this
##' function updates in order to define and fit the submodel.  This
##' should work for most standard linear-type models in common R packages. 
##'
##' @inheritParams fic
##'
##' @return List of all fitted submodel objects 
##'
##' @export
fit_submodels <- function(wide, inds){
    nmod <- nrow(inds)
    sub <- vector(nmod, mode="list")
    XZ <- model.matrix(wide)
    for (i in 1:nmod){
        XZi <- XZ[,which(inds[i,]==1),drop=FALSE]
        call <- wide$call
        call$formula <- update(as.formula(call$formula), . ~ XZi - 1)
        sub[[i]] <- eval(call)
    }
    sub
}
