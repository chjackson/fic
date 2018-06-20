##' Focused Information Criterion for parametric survival models
##'
##' Focused information criterion for parametric survival models fitted with the \pkg{survival} package.  
##'
##' @param wide Object returned by the \code{\link{survreg}} function, containing the wide model.  
##'
##' @param sub Object of class \code{\link{survreg}} containing the submodel to be assessed.   Optional. Only required if you want the estimate of the focus
##' function under the submodel to be included in the results. 
##' 
##' @details Any situation where all models being compared are special cases of a single "wide" model are supported.  Examples include covariate selection, selection between models for the baseline hazard/survival with different levels of flexibility (e.g. comparing exponential and Weibull).  Some of these are illustrated in the "survival" vignette.
##'
##' @inheritParams fic
##' 
##' @export
fic.survreg <- function(wide, inds, inds0=NULL, gamma0=0, focus=NULL, focus_deriv=NULL, X=NULL, sub=NULL, ...){
    coef_fn = function(x){        
        x$icoef
    }
    flexsurvreg_fns <- list(
        coef = function(x){x$icoef},
        nobs = function(x)nrow(model.frame(x))
    )
    fic.default(wide=wide, inds=inds, inds0=inds0, gamma0=gamma0, 
        focus=focus, focus_deriv=focus_deriv, X=X, sub=sub,
        fns = flexsurvreg_fns, ...)
}
