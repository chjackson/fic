##' Focused information criteria for parametric survival models
##'
##' Focused information criteria for parametric survival models fitted with the \pkg{survival} package.  
##'
##' @param wide Object returned by the \code{\link{survreg}} function, containing the wide model.  
##'
##' @param sub List of fitted model objects of class \code{\link{survreg}} containing the submodels to be assessed.   Optional. Only required if you want the estimate of the focus
##' function under the submodels to be included in the results. 
##' 
##' @details Any situation where all models being compared are special cases of a single "wide" model are supported.  Examples include covariate selection, selection between models for the baseline hazard/survival with different levels of flexibility (e.g. comparing exponential and Weibull). An example of the latter is in the \pkg{fic} package vignette "Examples of focused model comparison: parametric survival models".
##'
##' Parameters \code{par} of the focus function should be on the scale reported by the \code{icoef} component of the results of \code{survreg}, that is, with any positive-valued parameters log transformed.   
##'
##' @inheritParams fic
##'
##' @examples
##'
##' library(survival)
##' 
##' ## Fit exponential and Weibull models and plot fitted survival curves
##' ex <-  survreg(Surv(futime, fustat) ~ 1, data=ovarian, dist="exponential")
##' we <-  survreg(Surv(futime, fustat) ~ 1, data=ovarian, dist="weibull")
##' 
##' ## Plot fitted survival curves, highlighting 1 year survival 
##' plot(survfit(Surv(futime, fustat) ~ 1, data=ovarian))
##' t <- seq(0, 1200)
##' lines(t, pweibull(q=t, shape=exp(we$icoef[2]),
##'                   scale=exp(we$icoef[1]), lower.tail=FALSE))
##' lines(t, pexp(q=t, rate=1/exp(ex$icoef[1]), lower.tail=FALSE), lty=2)
##' abline(v=365, col="gray")
##' 
##' ## Focused model comparison for focus of 1-year survival probability
##' indmat <- rbind(exp    = c(1,0),
##'                 weib   = c(1,1))
##' surv1yr <- function(par){
##'     pweibull(q=365, shape=exp(par[2]), scale=exp(par[1]), lower.tail=FALSE)
##' }
##' fic(we, inds=indmat, focus=surv1yr, sub=list(ex, we))
##' 
##' ## Exponential model has lower expected error, given such a small dataset 
##' 
##' @export
fic.survreg <- function(wide, inds, inds0=NULL, gamma0=0, focus=NULL, focus_deriv=NULL, wt=NULL, sub=NULL, B=0, loss=loss_mse, ...){
    survreg_fns <- list(
        coef = function(x){x$icoef},
        nobs = function(x)nrow(model.frame(x))
    )
    fic.default(wide=wide, inds=inds, inds0=inds0, gamma0=gamma0, 
        focus=focus, focus_deriv=focus_deriv, wt=wt, sub=sub, B=B, FIC=TRUE, loss=loss, 
        fns = survreg_fns, ...)
}
