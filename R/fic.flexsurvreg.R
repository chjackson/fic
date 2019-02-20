##' Focused Information Criterion for flexible parametric survival models
##'
##' Focused information criterion for parametric survival models fitted with the \pkg{flexsurv} package.  
##'
##' @param wide Object of class \code{"flexsurvreg"} containing the wide model.  
##' These are returned, for example, by the \code{\link[flexsurv]{flexsurvreg}} and 
##' the \code{\link[flexsurv]{flexsurvspline}} functions.
##'
##' @param sub Object of class \code{\link[flexsurv]{flexsurvreg}} containing the submodel to be assessed.   Optional. Only required if you want the estimate of the focus
##' function under the submodel to be included in the results. 
##' 
##' @details Any situation where all models being compared are special cases of a single "wide" model are supported.  Examples include covariate selection, selection between models for the baseline hazard/survival with different levels of flexibility (e.g. comparing exponential, Weibull and generalized gamma).  Some of these are illustrated in the "survival" vignette.
##'
##' The choice between \code{\link[flexsurv]{flexsurvspline}} models with different numbers of knots is not supported, unless perhaps if the knot locations are defined manually so that models are nested within each other, but this has not been investigated.
##' 
##' @inheritParams fic
##' 
##' @export
fic.flexsurvreg <- function(wide, inds, inds0=NULL, gamma0=0, focus=NULL, focus_deriv=NULL, X=NULL, Xwt=NULL, sub=NULL, B=0, loss=loss_mse, ...){
    coef_fn = function(x){        
        if (!is.null(x$fixedpars))
            x$res.t[,"est"][-x$fixedpars]
        else x$res.t[,"est"]
    }
    flexsurvreg_fns <- list(
        coef = coef_fn,
        nobs = function(x)nrow(model.frame(x))
    )
    fic.default(wide=wide, inds=inds, inds0=inds0, gamma0=gamma0, 
                focus=focus, focus_deriv=focus_deriv, X=X, Xwt=Xwt, sub=sub,
                B=B, FIC=TRUE, loss=loss, 
                fns = flexsurvreg_fns, ...)
}
