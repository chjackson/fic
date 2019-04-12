##' Focused information criteria for linear models
##'
##' Focused information criteria for linear models fitted with \code{\link{lm}}.
##' Typically used to compare models with different covariates (more generally,
##' different linear terms).
##'
##' @param sub If \code{"auto"} (the default) then the submodels are fitted automatically within this function.   If \code{NULL} they are not fitted, and focus estimates are not returned with the results. 
##'
##' The model parameters include the intercept, followed by the coefficients of any covariates.   The standard deviation is excluded.
##'
##' Only covariate selection problems are supported in this function.  To compare between models with a fixed and unknown standard deviation, hand-written maximum likelihood estimation routines would be needed, along the lines described in the "skew-normal models" vignette. 
##'
##' The focus can depend on the standard deviation.  The focus function should then have an argument \code{sigma}. 
##'
##' See the vignette "Using the fic R package for focused model comparison: linear regression" for some examples. 
##' 
##' @inheritParams fic
##'
##' @examples
##'
##' ## Covariate selection in Motor Trend cars data
##' ## See the "fic" package vignette on linear models for more details
##'
##' wide.lm <- lm(mpg ~ am + wt + qsec + disp + hp, data=mtcars)
##'
##' ## Select between all submodels 
##' ncovs_wide <- length(coef(wide.lm)) - 1
##' inds0 <- c(1, rep(0, ncovs_wide))
##' inds <- all_inds(wide.lm, inds0)
##'
##' ## Two focuses: mean MPG for automatic and manual transmission,
##' ## given mean values of the other covariates 
##' cmeans <- colMeans(model.frame(wide.lm)[,c("wt","qsec","disp","hp")])
##' X <- rbind(
##'   "auto"   = c(intercept=1, am=0, cmeans),
##'   "manual" = c(intercept=1, am=1, cmeans)
##' )
##' ficres <- fic(wide.lm, inds=inds, focus=mean_normal, X=X)
##' summary(ficres)
##' ggplot_fic(ficres)
##'
##'
##' @export
fic.lm <- function(wide, inds, inds0=NULL, gamma0=0, focus=NULL, focus_deriv=NULL,
                          wt=NULL, sub="auto", B=0, loss=loss_mse, ...){
    if (!inherits(wide, "lm")) stop("\"wide\" must be an object of class \"lm\"")
    fns <- list(
        aux = function(x)list(sigma=summary(x)$sigma)
    )
    fic.default(wide=wide, inds=inds, inds0=inds0, gamma0=gamma0, 
                focus=focus, focus_deriv=focus_deriv, wt=wt, sub=sub,
                fns=fns, B=B, FIC=TRUE, loss=loss, ...)
}

