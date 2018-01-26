##' Focused Information Criterion for flexible parametric survival models
##'
##' Focused information criterion for parametric survival models fitted with the \pkg{flexsurv} package.  
##'
##' @param wide Object of class \code{\link{"flexsurvreg"}} containing the wide model.  
##' These are returned, for example, by the \code{\link{flexsurvreg}} and 
##' the \code{\link{flexsurvspline}} functions.
##'
##' @param sub Object of class \code{\link{flexsurvreg}} containing the submodel to be assessed.  
##'
##' @param inds TODO better documentation for what indices correspond to what parameters, and what model selection problems are supported.
##' Do at least covariate selection, selection between gamma, generalised gamma etc.
##' Clarify choice between splines with different numbers of knots not supported, unless defined specially with nested knot locations. 
##' 
##' @inheritParams fic
##' 
##' @export
fic.flexsurvreg <- function(wide, sub=NULL, inds, inds0, gamma0=0, focus=NULL, focus_deriv=NULL, X=NULL, ...){
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


