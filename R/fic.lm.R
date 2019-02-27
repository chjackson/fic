##' Focused Information Criterion for linear models
##'
##' Focused information criterion for linear models fitted with \code{\link{lm}}.
##' Used to compare models with different covariates (more generally,
##' different linear terms).
##'
##' @param sub If \code{"auto"} (the default) then the submodels are fitted automatically within this function.   If \code{NULL} they are not fitted, and focus estimates are not returned with the results. 
##'
##' TODO put mtcars in rd example 
##' 
##' @inheritParams fic
##'
##' @export
fic.lm <- function(wide, inds, inds0=NULL, gamma0=0, focus=NULL, focus_deriv=NULL,
                          X=NULL, Xwt=NULL, sub="auto", B=0, loss=loss_mse, ...){
    if (!inherits(wide, "lm")) stop("\"wide\" must be an object of class \"lm\"")
    fic.default(wide=wide, inds=inds, inds0=inds0, gamma0=gamma0, 
                focus=focus, focus_deriv=focus_deriv, X=X, Xwt=Xwt, sub=sub,
                B=B, FIC=TRUE, loss=loss, ...)
}

