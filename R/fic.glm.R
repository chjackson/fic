##' Focused Information Criterion for generalized linear models
##'
##' Focused information criterion for generalized linear models fitted with glm().
##' This can be used, for example, to compare models with different covariates.
##'
##' @inheritParams fic
##'
##' @export
fic.glm <- function(wide, inds, inds0=NULL, gamma0=0, focus=NULL, focus_deriv=NULL,
                          X=NULL, Xwt=NULL, sub=NULL, B=0, loss=loss_mse, ...){
    if (!inherits(wide, "glm")) stop("\"wide\" must be an object of class \"glm\"")
    fic.default(wide=wide, inds=inds, inds0=inds0, gamma0=gamma0, 
                focus=focus, focus_deriv=focus_deriv, X=X, Xwt=Xwt, sub=sub, B=B, loss=loss, ...)
}

