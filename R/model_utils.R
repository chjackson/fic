##' Form `fic` model indicator argument in presence of factors 
##'
##' Given a model indicator \code{inds} identifying terms in a regression
##' model, convert this to the format needed for \code{\link{fic}} by converting
##' indicators for regression terms to indicators for inclusion of
##' parameters.  Only required if there are factors.
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
##' @examples
##'
##' # Five terms in this model: intercept and four covariates,
##' # but the covariate "ftv" is a factor with 3 levels,
##' # so there are six parameters
##' bwt.glm <- glm(low ~ lwtkg + age + smoke + ftv, data=birthwt, family="binomial")
##'
##' ## Convert indicator for terms to indicator for parameters 
##' inds <- rbind(c(1,1,1,0,0),
##'               c(1,1,1,1,1))
##' expand_inds(inds, bwt.glm)
##' 
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
            stop(sprintf("`inds` has %s columns, but %s terms in model%s", ncol(inds), nterms, icpt_msg))
        }
    }
    inds[,ass,drop=FALSE]
}



##' Form indicator matrix describing all submodels of a general linear wide model
##'
##' Form indicator matrix describing all possible submodels of a general linear wide model, where the submodels are defined by selected covariates.
##'
##' @param wide A fitted model of standard R format, such that \code{terms(wide)} returns information about the terms of the model formula.  Models outside standard R packages may not support this.
##'
##' @param inds0 Narrow model indicators, in format described in \code{\link{fic}}.
##'
##' @param auxpars Names of parameters in the wide model other than the covariate effects being selected from.  By default, for linear and generalised linear models this is \code{c("(Intercept)")}, and for Cox regression this is omitted.
##'
##' @param ... Other arguments. Currently unused. 
##'
##' @return A matrix in the format required by the \code{inds} argument of \code{fic()}, representing all possible submodels of the wide model.
##'
##' The number of rows is the number of models, and the number of columns is the number of parameters in the wide model.   The \eqn{r,s} entry of the matrix is a 1 if the \eqn{r}th submodel includes parameter \eqn{s}, and 0 otherwise.
##'
##' If a factor is included (excluded) from the submodel, then all corresponding parameters are included (excluded).
##'
##' @examples
##' bwt.glm <- glm(low ~ lwtkg + age + smoke, data=birthwt, family="binomial")
##' all_inds(bwt.glm, inds0=c(1,0,0,0))
##'
##' # note no intercept term in Cox models, so inds0 has two elements here
##' library(survival)
##' wide <- coxph(Surv(years, death==1) ~ sex + thick_centred, data=melanoma)
##' all_inds(wide, inds0=c(0,0))
##' 
##' @rdname all_inds 
##' 
##' @export
all_inds.default <- function(wide, inds0=NULL, auxpars=NULL, ...){
    tt <- terms(wide)
    labs <- attr(tt, "term.labels")
    labs <- c(auxpars,labs)
    inds_short <- do.call(expand.grid, Map(function(x)c(0,1), labs))
    rownames(inds_short) <- apply(inds_short, 1, paste, collapse="")
    asn <- attr(model.matrix(wide), "assign")
    pick <- match(asn, unique(asn))
    combs <- inds_short[,pick,drop=FALSE]
    if(length(inds0) != ncol(combs))
        stop(sprintf("inds0 of length %s, but %s parameters in wide model",
                     length(inds0),ncol(combs)))
    ininds0 <- if (!is.null(inds0)) apply(combs, 1, function(x) !any(x==0 & inds0==1)) else rep(TRUE, nrow(combs))
    combs <- combs[ininds0,]
    attr(combs, "termnames") <- labs
    combs
}

##' @rdname all_inds 
##' @export
all_inds <- function(wide,inds0,...) UseMethod("all_inds")

## todo put in model specific files

##' @rdname all_inds 
##' @export
all_inds.lm <- function(wide, inds0, ...) {
    all_inds.default(wide=wide, inds0=inds0, auxpars=c("(Intercept)"))
}

##' @rdname all_inds 
##' @export
all_inds.glm <- function(wide, inds0, ...) {
    all_inds.default(wide=wide, inds0=inds0, auxpars=c("(Intercept)"))
}

##' @rdname all_inds 
##' @export
all_inds.coxph <- function(wide, inds0, ...) {
    all_inds.default(wide=wide, inds0=inds0, auxpars=NULL)
}


##' Convert data frame of covariate values to a design matrix
##'
##' @param newdata Data frame where each row is a vector of covariate values defining an alternative focus quantity.
##'
##' @param wide Wide model which includes these covariates.
##'
##' @param intercept Include an intercept as the first column.
##'
##' @return "Design" matrix of covariate values defining alternative focuses, with factors expanded to their contrasts.  This is in the form required by the \code{X} argument of \code{\link{fic}}, with one row per alternative focus. The columns correspond to coefficients in a linear-type model.  For the built-in focus functions such as \code{\link{mean_normal}} and \code{\link{prob_logistic}}, these coefficients include an intercept, but user-written focuses may be written in such a way as not to require an intercept (as in the example in the "skew normal" vignette). 
##'
##' @details Numeric values can be supplied for factor levels that are character strings denoting numbers (like \code{"1"} or \code{"2"}).
##'
##' @examples
##' bwt.glm <- glm(low ~ lwtkg + age + smoke + ftv, data=birthwt, family="binomial")
##' newdata <- data.frame(lwtkg=1, age=60, smoke=0, ftv="2+")
##' newdata_to_X(newdata, bwt.glm)
##' 
##' ## See the Cox regression section of the main package vignette for another example.
##'
##' @export
newdata_to_X <- function(newdata, wide, intercept=TRUE){
    tt <- terms(wide)
    faclevs <- .getXlevels(tt, model.frame(wide)) # list of levels of all factors
    fn <- names(faclevs)
    newdata[fn] <- lapply(newdata[fn], as.character) ## convert numerics supplied for factor levels to characters
    X <- model.matrix(delete.response(tt), data=newdata, xlev=faclevs)
    if (!intercept) X <- X[,-1,drop=FALSE]
    X
}


##' Fit submodels of a general linear wide model, defined by a matrix of indicators
##' for inclusion of covariates
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
##' @return List of all fitted submodel objects.
##'
##' @rdname fit_submodels
##'
##' @examples
##' bwt.glm <- glm(low ~ lwtkg + age + smoke,
##'                data=birthwt, family="binomial")
##' inds <- rbind(c(1,1,1,0), c(1,1,0,0))
##' fit_submodels(bwt.glm, inds=inds)
##' 
##' @export
fit_submodels <- function(wide,inds,...) UseMethod("fit_submodels")

##' @export
fit_submodels.default <- function(wide, inds, ...){
    nmod <- nrow(inds)
    sub <- vector(nmod, mode="list")
    names(sub) <- rownames(inds)
    XZ <- model.matrix(wide)
    for (i in 1:nmod){
        XZi <- XZ[,which(inds[i,]==1),drop=FALSE]
        call <- wide$call
        call$formula <- update(as.formula(call$formula), . ~ XZi - 1)
        sub[[i]] <- eval(call)
    }
    sub
}


### Can't seem to just do fit_submodels and then extract basehaz
### externally, since then basehaz seems to see the wrong XZi, for some
### strange reason related to environments.

##' @export
fit_submodels.coxph <- function(wide, inds, ...){
    nmod <- nrow(inds)
    sub <- vector(nmod, mode="list")
    names(sub) <- rownames(inds)
    XZ <- model.matrix(wide)
    for (i in 1:nmod){
        XZi <- XZ[,which(inds[i,]==1),drop=FALSE]
        call <- wide$call
        call$formula <- update(as.formula(call$formula), . ~ XZi - 1)
        sub[[i]] <- eval(call)
        attr(sub[[i]], "basehaz") <- basehaz(sub[[i]], centered=FALSE)
    }
    sub
}
