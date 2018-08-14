##' Focused Information Criterion for linear models
##'
##' Focused information criterion for linear models fitted with lm().
##' Used to compare models with different covariates (more generally,
##' different linear terms).
##'
##' (something about sub argument)
##'
##' TODO put mtcars in rd example 
##' 
##' @inheritParams fic
##'
##' @export
fic.lm <- function(wide, inds, inds0=NULL, gamma0=0, focus=NULL, focus_deriv=NULL,
                          X=NULL, Xwt=NULL, sub="auto", B=0, loss=loss_mse, ...){
    if (!inherits(wide, "lm")) stop("\"wide\" must be an object of class \"lm\"")
    fic.default(wide=wide, inds=inds, inds0=inds0, gamma0=gamma0, 
                focus=focus, focus_deriv=focus_deriv, X=X, Xwt=Xwt, sub=sub, B=B, loss=loss, ...)
}

