##' Focused information criteria for flexible parametric survival models
##'
##' Focused information criteria for parametric survival models fitted with the \pkg{flexsurv} package.  
##'
##' @param wide Object of class \code{"flexsurvreg"} containing the wide model.  
##' These are returned, for example, by the \code{\link[flexsurv]{flexsurvreg}} and 
##' the \code{\link[flexsurv]{flexsurvspline}} functions.
##'
##' @param sub List of objects of class \code{\link[flexsurv]{flexsurvreg}} containing the submodels to be assessed. Optional. Only required if you want the estimate of the focus
##' function under the submodels to be included in the results. 
##' 
##' @details Any situation where all models being compared are special cases of a single "wide" model are supported.  Examples include covariate selection, selection between models for the baseline hazard/survival with different levels of flexibility (e.g. comparing exponential, Weibull and generalized gamma).  Some of these are illustrated in the \pkg{fic} package vignette "Examples of focused model comparison: parametric survival models".
##'
##' The choice between \code{\link[flexsurv]{flexsurvspline}} models with different numbers of knots is not supported, unless perhaps if the knot locations are defined manually so that models are nested within each other, but this has not been investigated.
##' 
##' @inheritParams fic
##' 
##' @export
##'
##' @examples
##'
##' ## Simulated example from the "fic" package vignette on
##' ## parametric survival modelling.
##' ## See this vignette for more details and more examples.
##' 
##' set.seed(1)
##'
##' if (requireNamespace("flexsurv", quietly=TRUE)){
##'
##' ## Simulate from an exponential 
##' y <- rexp(50); cen <- rep(1,50)
##'
##' ## Fit wide generalized gamma, and compare
##' ## exponential, weibull and generalized gamma models
##' indmat <- rbind(exp    = c(1,0,0),
##'                 weib   = c(1,1,0),
##'                 ggamma = c(1,1,1))
##' gge <- flexsurv::flexsurvreg(survival::Surv(y, cen) ~ 1, dist="gengamma")
##'
##' ## Focus is restricted mean survival over 8 time units 
##' focus <- function(par){
##'    flexsurv::rmst_gengamma(8, par[1], exp(par[2]), par[3])
##' }
##'
##' ## Weibull model actually has lowest FIC and RMSE even though it's
##' ## not true: extra variability is deemed worth alleviating the
##' ## risk of bias.
##' 
##' fic(gge, inds=indmat, gamma0=c(0,1), focus=focus)
##'
##' }


fic.flexsurvreg <- function(wide, inds, inds0=NULL, gamma0=0, focus=NULL, focus_deriv=NULL, wt=NULL, sub=NULL, B=0, loss=loss_mse, ...){
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
                focus=focus, focus_deriv=focus_deriv, wt=wt, sub=sub,
                B=B, FIC=TRUE, loss=loss, 
                fns = flexsurvreg_fns, ...)
}
