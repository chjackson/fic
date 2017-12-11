##' Focused Information Criterion for multi-state models for panel data
##'
##' Focused information criterion for multi-state models fitted with msm().  This is not yet a S3 method, but will probably be made into one
##'
##' @param object Object returned by \code{\link{msm}} containing the fitted model.
##'
##' @param inds TODO better documentation for what indices correspond to what parameters, and what model selection problems are supported.  Just covariate selection?  what about constraints?
##' 
##' @inheritParams fic
##' 
##' @export
fic.msm <- function(object, inds, pp, focus=NULL, focus_deriv=NULL){
    ests <- object$estimates
    n <- attr(model.frame(object), "ntrans")
    J <- object$opt$hessian / n # TODO or use expected information, will this be more accurate?
    fic(ests=ests, J=J, inds=inds, pp=pp, n=n,
        focus=focus, focus_deriv=focus_deriv)
}


