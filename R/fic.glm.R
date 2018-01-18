##' Focused Information Criterion for generalized linear models
##'
##' Focused information criterion for generalized linear models fitted with glm().
##'
##' @param wide Object returned by \code{\link{glm}} containing the wide model.
##'
##' @param sub Object returned by \code{\link{glm}} containing the submodel to be assessed.  Optional. 
##' Only required if you want the estimate of the focus function under the submodel to be included in the results. 
##' 
##' @inheritParams fic
##' 
##' @export
fic.glm <- function(wide, sub=NULL, inds, inds0, focus=NULL, focus_deriv=NULL, ...){
    par <- coef(wide)
    n <- nobs(wide)
    J <- solve(vcov(wide)) / n
    fic(par=par, J=J, inds=inds, inds0=inds0, n=n,
              focus=focus, focus_deriv=focus_deriv, 
              parsub=coef(sub), ...)
}

