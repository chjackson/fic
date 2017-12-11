##' Focused Information Criterion for generalized linear models
##'
##' Focused information criterion for generalized linear models fitted with glm().  This is not yet a S3 method, but will probably be made into one
##'
##' @param object Object returned by \code{\link{glm}} containing the fitted model.
##' 
##' @inheritParams fic
##' 
##' @export
fic.glm <- function(object, inds, pp, focus=NULL, focus_deriv=NULL, ...){
    ests <- coef(object)
    n <- nobs(object)
    J <- solve(vcov(object)) / n
    fl <- get_focus(focus, focus_deriv, ests, ...)
    fic(ests=ests, J=J, inds=inds, pp=pp, n=n,
        focus=fl$focus, focus_deriv=fl$focus_deriv)
}

