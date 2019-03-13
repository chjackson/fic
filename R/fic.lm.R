##' Focused information criteria for linear models
##'
##' Focused information criteria for linear models fitted with \code{\link{lm}}.
##' Typically used to compare models with different covariates (more generally,
##' different linear terms).
##'
##' @param sub If \code{"auto"} (the default) then the submodels are fitted automatically within this function.   If \code{NULL} they are not fitted, and focus estimates are not returned with the results. 
##'
##' The model parameters include the intercept, followed by the coefficients of any covariates.   The standard deviation is excluded.
##'
##' Only covariate selection problems are supported in this function.  To compare between models with a fixed and unknown standard deviation, hand-written maximum likelihood estimation routines would be needed, along the lines described in the "skew-normal models" vignette. 
##'
##' The focus can depend on the standard deviation.  The focus function should then have an argument \code{sigma}. 
##'
##' See the vignette "Using the fic R package for focused model comparison: linear regression" for some examples. 
##' 
##' @inheritParams fic
##'
##' @export
fic.lm <- function(wide, inds, inds0=NULL, gamma0=0, focus=NULL, focus_deriv=NULL,
                          X=NULL, Xwt=NULL, sub="auto", B=0, loss=loss_mse, ...){
    if (!inherits(wide, "lm")) stop("\"wide\" must be an object of class \"lm\"")
    fns <- list(
        aux = function(x)list(sigma=summary(x)$sigma)
    )
    fic.default(wide=wide, inds=inds, inds0=inds0, gamma0=gamma0, 
                focus=focus, focus_deriv=focus_deriv, X=X, Xwt=Xwt, sub=sub,
                fns=fns, B=B, FIC=TRUE, loss=loss, ...)
}

