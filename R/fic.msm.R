##' Focused information criteria for multi-state models for panel data
##'
##' Focused information criteria for multi-state models fitted with \code{\link[msm]{msm}} from the \pkg{msm} package. 
##'
##' @param wide Object returned by \code{\link[msm]{msm}} containing the wide model.
##'
##' @param sub List of objects returned by \code{\link[msm]{msm}} containing the submodels to be assessed.  Optional. Only required if you want the estimate of the focus
##' function under the submodel to be included in the results. 
##'
##' @details This might be used for covariate selection, or comparing models with different constraints on the covariate effects or intensities.   An example is given in the \pkg{fic} package vignette "Examples of focused model comparison: multi-state models".  Note in particular in this example how the parameters are ordered in the \code{inds} argument, and how the various \pkg{msm} output functions can be used as focuses. 
##' 
##' @inheritParams fic
##'
##' @examples
##'
##' ## Covariate selection in psoriatic arthritis model.
##' ## See the "fic" package vignette on multi-state models for
##' ## more details and examples. 
##'
##' if (requireNamespace("msm",quietly=TRUE)){
##'
##' Qind <- rbind(c(0, 1, 0, 0),
##'               c(0, 0, 1, 0),
##'               c(0, 0, 0, 1),
##'               c(0, 0, 0, 0))
##' psor.wide.msm <- msm::msm(state ~ months, subject=ptnum, data=msm::psor, 
##'                      qmatrix = Qind,  gen.inits=TRUE,
##'                      covariates = ~ollwsdrt+hieffusn)
##' inds <- rbind(
##'     c(1,1,1,0,0,0,0,0,0),
##'     c(1,1,1,0,0,0,0,0,1),
##'     c(1,1,1,0,0,0,0,1,1),
##'     c(1,1,1,0,0,0,1,1,1),
##'     c(1,1,1,0,0,1,1,1,1),
##'     c(1,1,1,0,1,1,1,1,1),
##'     c(1,1,1,1,1,1,1,1,1)
##' )
##' focus_tlos <- function(par){
##'     x.new <- msm::updatepars.msm(psor.wide.msm, par)
##'     msm::totlos.msm(x.new, covariates=0, tot=10)["State 4"]
##' }
##' fres <- fic(wide=psor.wide.msm, inds=inds, focus=focus_tlos)
##' fres
##' 
##' }
##' 
##' @export
fic.msm <- function(wide, inds, inds0=NULL, gamma0=0, focus=NULL, focus_deriv=NULL, wt=NULL, sub=NULL, B=0, loss=loss_mse, ...){
    msm_fns <- list(coef = function(x)x$estimates,
                    nobs = function(x){
                       mf <- model.frame(x)
                       nrow(mf) - attr(mf, "npts") # attr(mf, "ntrans") in recent versions
                    },
                    vcov = function(x)x$covmat)
    fic.default(wide=wide, inds=inds, inds0=inds0, gamma0=gamma0, 
                focus=focus, focus_deriv=focus_deriv, wt=wt, sub=sub,
                B=B, FIC=TRUE, loss=loss, 
                fns = msm_fns, ...)
}
