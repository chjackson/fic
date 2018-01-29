##' Focused Information Criterion for multi-state models for panel data
##'
##' Focused information criterion for multi-state models fitted with msm().
##'
##' @param wide Object returned by \code{\link{msm}} containing the wide model.
##'
##' @param sub Object returned by \code{\link{msm}} containing the submodel to be assessed.  Optional. Only required if you want the estimate of the focus
##' function under the submodel to be included in the results. 
##'
##' @param inds TODO better documentation for what indices correspond to what parameters, and what model selection problems are supported.  Just covariate selection?  what about constraints?
##' 
##' @inheritParams fic
##' 
##' @export
fic.msm <- function(wide, inds, inds0, gamma0=0, focus=NULL, focus_deriv=NULL, X=NULL, sub=NULL, ...){
    msm_fns <- list(coef = function(x)x$estimates,
                    nobs = function(x)attr(model.frame(x), "ntrans"),
                    vcov = function(x)x$covmat)
    fic.default(wide=wide, inds=inds, inds0=inds0, gamma0=gamma0, 
        focus=focus, focus_deriv=focus_deriv, X=X, sub=sub,
        fns = msm_fns, ...)
}
