##' Focused information criteria for multi-state models for panel data
##'
##' Focused information criteria for multi-state models fitted with \code{\link[msm]{msm}} from the \pkg{msm} package. 
##'
##' @param wide Object returned by \code{\link[msm]{msm}} containing the wide model.
##'
##' @param sub List of objects returned by \code{\link[msm]{msm}} containing the submodels to be assessed.  Optional. Only required if you want the estimate of the focus
##' function under the submodel to be included in the results. 
##'
##' @details This might be used for covariate selection, or comparing models with different constraints on the covariate effects or intensities.   An example is given in the \pkg{fic} package vignette "Examples of focused model comparison: multi-state models".  Note in particular in this example how the parameters are ordered in the \code{inds} argument, and how the various \pkg{msm} output functions can be used as focuses. 
##' 
##' @inheritParams fic
##' 
##' @export
fic.msm <- function(wide, inds, inds0=NULL, gamma0=0, focus=NULL, focus_deriv=NULL, wt=NULL, sub=NULL, B=0, loss=loss_mse, ...){
    msm_fns <- list(coef = function(x)x$estimates,
                    nobs = function(x){
                       mf <- model.frame(x)
                       nrow(mf) - attr(mf, "npts") # attr(mf, "ntrans") in recent versions
                    },
                    vcov = function(x)x$covmat)
    fic.default(wide=wide, inds=inds, inds0=inds0, gamma0=gamma0, 
                focus=focus, focus_deriv=focus_deriv, wt=wt, sub=sub,
                B=B, FIC=TRUE, loss=loss, 
                fns = msm_fns, ...)
}
