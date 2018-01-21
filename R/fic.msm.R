##' Focused Information Criterion for multi-state models for panel data
##'
##' Focused information criterion for multi-state models fitted with msm().  This is not yet a S3 method, but will probably be made into one
##'
##' @param wide Object returned by \code{\link{msm}} containing the wide model.
##'
##' @param sub Object returned by \code{\link{msm}} containing the submodel to be assessed.  
##'
##' @param inds TODO better documentation for what indices correspond to what parameters, and what model selection problems are supported.  Just covariate selection?  what about constraints?
##' 
##' @inheritParams fic
##' 
##' @export
fic.msm <- function(wide, sub=NULL, inds, inds0, focus=NULL, focus_deriv=NULL, X=NULL, ...){
    par <- wide$estimates
    n <- attr(model.frame(wide), "ntrans")
    J <- wide$opt$hessian / n # or use expected information, will this be more accurate?
    fic(par=par, J=J, inds=inds, inds0=inds0, n=n,
               focus=focus, focus_deriv=focus_deriv, 
                parsub=sub$estimates, X=X, ...)
}


