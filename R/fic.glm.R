##' Focused Information Criterion for generalized linear models
##'
##' Focused information criterion for generalized linear models fitted with \code{\link{glm}}.
##' Used to compare models with different covariates (more generally,
##' different linear terms).
##'
##' @param sub If \code{"auto"} (the default) then the submodels are fitted automatically within this function.   If \code{NULL} they are not fitted, and focus estimates are not returned with the results. 
##'
##' @inheritParams fic
##'
##' @export
fic.glm <- function(wide, inds, inds0=NULL, gamma0=0, focus=NULL, focus_deriv=NULL,
                          X=NULL, Xwt=NULL, sub="auto", B=0, loss=loss_mse, ...){
    if (!inherits(wide, "glm")) stop("\"wide\" must be an object of class \"glm\"")
    fic.default(wide=wide, inds=inds, inds0=inds0, gamma0=gamma0, 
                focus=focus, focus_deriv=focus_deriv, X=X, Xwt=Xwt,
                sub=sub, B=B, FIC=TRUE, loss=loss, ...)
}

