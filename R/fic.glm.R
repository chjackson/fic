##' Focused information criteria for generalized linear models
##'
##' Focused information criteria for generalized linear models fitted with \code{\link{glm}}.
##' Used to compare models with different covariates (more generally,
##' different linear terms).
##'
##' @param sub If \code{"auto"} (the default) then the submodels are fitted automatically within this function.   If \code{NULL} they are not fitted, and focus estimates are not returned with the results. 
##'
##' @inheritParams fic
##'
##' @details 
##' The model parameters include the intercept, followed by the coefficients of any covariates.  Any "dispersion" parameter is excluded from the parameters indicated by \code{inds} and \code{inds0}.   
##'
##' Only covariate selection problems are supported in this function.
##' To compare between models with a fixed and unknown dispersion,
##' \code{\link{glm}} would have to be replaced by maximum likelihood
##' estimation routines written by hand, along the lines described in
##' the "skew-normal models" vignette.
##'
##' The focus function can however depend on the value of the
##' dispersion parameter.  The focus function should then have
##' an argument called \code{dispersion}.  See the example of a gamma
##' GLM below, where the focus is the mean outcome for some covariate
##' value.
##'
##' Examples of covariate selection in logistic regression are given in the main package vignette.
##'
##' @examples
##'
##' # Gamma regression with one binary covariate
##'
##' # Simulated data 
##' set.seed(1)
##' simx <- rbinom(1000, 1, 0.5)
##' mean0 <- 1.1
##' simy <- rgamma(1000, shape=2, scale=exp(log(mean0) + 0.2*simx))
##' 
##' mod <- glm(simy ~ simx, family=Gamma(link="log"))
##'
##' # Check the parameter estimates are close to true
##' # values used for simulation
##' (shape <- 1 / summary(mod)$dispersion)
##' coef(mod)[2] # log mean ratio associated with covariate. true value 0.2
##' exp(coef(mod)[1] - log(shape))  # mean with x=0, true value 1.1
##' exp(coef(mod)[1] + coef(mod)[2] - log(shape))  # mean with x=1
##'
##' focus_mean <- function(par, X, dispersion){  
##'   exp(X %*% par - log(1/dispersion))
##' }
##'
##' X <- rbind("x0" = c(1,0), "x1" = c(1,1))
##' inds <- rbind("no_covariate"=c(1,0), "covariate"=c(1,1))
##'
##' fic(mod, inds=inds, focus=focus_mean, X=X)
##'
##' # The focus need not depend on X or the dispersion
##' focus_base_scale <- function(par, dispersion){  
##'   exp(par[1])
##' }
##' fic(mod, inds=inds, focus=focus_base_scale)
##'
##' # ...equivalently,
##' focus_base_scale <- function(par){  
##'   exp(par[1])
##' }
##' fic(mod, inds=inds, focus=focus_base_scale)
##'
##' 
##' 
##' @export
fic.glm <- function(wide, inds, inds0=NULL, gamma0=0, focus=NULL, focus_deriv=NULL,
                          wt=NULL, sub="auto", B=0, loss=loss_mse, ...){
    if (!inherits(wide, "glm")) stop("\"wide\" must be an object of class \"glm\"")
    fns <- list(
        aux = function(x)list(dispersion=summary(x)$dispersion)
    )
    fic.default(wide=wide, inds=inds, inds0=inds0, gamma0=gamma0, 
                focus=focus, focus_deriv=focus_deriv, wt=wt,
                sub=sub, fns=fns, B=B, FIC=TRUE, loss=loss, ...)
}

