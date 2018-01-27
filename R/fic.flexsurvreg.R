##' Focused Information Criterion for flexible parametric survival models
##'
##' Focused information criterion for parametric survival models fitted with the \pkg{flexsurv} package.  
##'
##' @param wide Object of class \code{\link{"flexsurvreg"}} containing the wide model.  
##' These are returned, for example, by the \code{\link{flexsurvreg}} and 
##' the \code{\link{flexsurvspline}} functions.
##'
##' @param sub Object of class \code{\link{flexsurvreg}} containing the submodel to be assessed.   Optional. Only required if you want the estimate of the focus
##' function under the submodel to be included in the results. 
##'
##' @param inds TODO better documentation for what indices correspond to what parameters, and what model selection problems are supported.
##' Do at least covariate selection, selection between gamma, generalised gamma etc.
##' 
##' The choice between \code{\link{flexsurvspline}} models with different numbers of knots is not supported, unless perhaps if the knot locations are defined manually so that models are nested within each other, but this has not been investigated.
##' 
##' @inheritParams fic
##' 
##' @export
fic.flexsurvreg <- function(wide, inds, inds0, gamma0=0, focus=NULL, focus_deriv=NULL, X=NULL, sub=NULL, ...){
    par <- coef(wide)
    n <- nrow(model.frame(wide)) # TODO add nobs method
    J <- solve(vcov(wide)) / n # or use expected information, will this be more accurate?
    res <- fic(par=par, J=J, inds=inds, inds0=inds0, gamma0=gamma0, n=n,
               focus=focus, focus_deriv=focus_deriv, 
                parsub=coef(sub), X=X, ...)
    if (!is.null(sub)){
      res <- cbind(res, AIC=AIC(sub), BIC=BIC(sub))
    }
    res
}


