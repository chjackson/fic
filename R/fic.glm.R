##' Focused Information Criterion for generalized linear models
##'
##' Focused information criterion for generalized linear models fitted with glm().  This is not yet a S3 method, but will probably be made into one
##'
##' @param object Object returned by \code{\link{glm}} containing the fitted model.
##'
##' @inheritParams fic
##' 
##' @export
fic.glm <- function(object, inds, pp, focus=NULL, focus.deriv=NULL){
    ests <- coef(object)
    n <- nobs(object)
    J <- solve(vcov(object)) / n
    fic(ests=ests, J=J, inds=inds, pp=pp, n=n,
        focus=focus, focus.deriv=focus.deriv)
}
