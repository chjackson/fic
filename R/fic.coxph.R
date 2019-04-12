##' Focused information criteria for Cox proportional hazard regression models
##'
##' Focused model comparison for Cox models fitted with \code{coxph} from the \code{survival} package.  Built-in focuses include the hazard ratio, survival and cumulative hazard. 
##' 
##' @aliases fic.coxph
##' 
##' @inheritParams fic.default
##'
##' @param X Covariate values to evaluate the focus at.  See \code{\link{fic}}, argument \code{...}, for the required format. 
##'
##' @param t times to evaluate the focus at.  Only relevant for survival and cumulative hazard focuses, as the hazard ratio is constant through time.
##'
##'
##' @param focus Three built-in focus quantities are supported:
##'
##' \code{"hr"} for the hazard ratio between individuals with covariate values of \code{X} and 0. 
##'
##' \code{"survival"} for the survival probability at time or times given in \code{t}.
##'
##' \code{"cumhaz"} for the cumulative hazard at time or times given in \code{t}.
##'
##' Alternatively, a list of three R functions can be supplied, with components named \code{"focus"}, \code{"deriv"} and \code{"dH"} respectively giving the focus, derivative with respect to the log hazard ratios, and derivative with respect to the cumulative hazards.   Each function should have arguments \code{par}, \code{H0}, \code{X} and \code{t}, giving the log hazard ratios, baseline cumulative hazard, covariate values and time points at which the focus function should be evaluated.  See the section "User-defined focuses" below for a little more information.
##'
##'
##' @param sub If \code{"auto"} (the default) then the submodels are fitted automatically within this function.   If \code{NULL} they are not fitted, and focus estimates are not returned with the results.
##'
##' @details Stratified Cox models are not currently supported.
##'
##' @examples
##'
##' ## Example of focused covariate selection in a Cox model
##' ## For more information see the main package vignette.
##' 
##' library(survival)
##' wide <- coxph(Surv(years, death==1) ~ sex + thick_centred + infilt +
##'               epith + ulcer + depth + age, data=melanoma)
##'
##' ## Define models to be selected from
##' ## Sex included in all models
##' ## Select between models including all combinations of other covariates
##' inds0 <- expand_inds(c(1,0,0,0,0,0,0), wide)
##' combs <- all_inds(wide, inds0)
##'
##' ## Covariate values defining focus 
##' newdata <- with(melanoma,
##'                 data.frame(sex = c("female","male"),
##'                            thick_centred = tapply(thick_centred, sex, mean),
##'                            infilt=4, epith=1, ulcer=1, depth=2,
##'                            age = tapply(age, sex, mean)))
##' X <- newdata_to_X(newdata, wide, intercept=FALSE)
##'
##' ## Focus is 5-year survival for these covariate values
##' ficall <- fic(wide, inds=combs, inds0=inds0, focus="survival", X=X, t=5)
##' ggplot_fic(ficall)
##' summary(ficall)
##' 
##' @rdname fic.coxph
##'
##' @import abind
##' @import tensor
##' @importFrom survival basehaz
##'
##' @section User-defined focuses:
##'
##' Each function should have four arguments:
##'
##' \code{par} Vector of estimated coefficients, the log hazard ratios in the Cox model.
##'
##' \code{H0}  Cumulative hazard estimate at a set of times, in the form of the output from \code{\link[survival]{basehaz}}.   The function \code{\link{get_H0}} can be used on this estimate to obtain the estimate at any other times by interpolation.
##'
##' \code{X}  Matrix of covariates, with \code{ncov} rows and \code{npar} columns, where \code{ncov} is the number of alternative covariate values definining alternative focuses we want to compare models for, and \code{npar} is the number of coefficients in the model.
##'
##' \code{t}  Vector of times defining alternative focus quantities (such as the survival)
##'
##' For examples, examine the source for the built-in functions
##'
##' \code{fic:::cox_hr,fic:::cox_hr_deriv,fic:::cox_hr_dH} for the hazard ratio between \code{X} and \code{0}
##'
##' \code{fic:::cox_cumhaz,fic:::cox_cumhaz_deriv,fic:::cox_cumhaz_dH} for the cumulative hazard
##'
##' \code{fic:::cox_survival,fic:::cox_survival_deriv,fic:::cox_survival_dH} for the survival.
##'
##' 
##' @export
fic.coxph <- function(wide, inds, inds0=NULL, gamma0=0,
                      focus, X=NULL, t=NULL, sub="auto", tidy=TRUE, ...)
{
    if (!inherits(wide, "coxph")) stop("\"wide\" must be an object of class \"coxph\"")
    if (!is.null(attr(wide$terms, "specials")$strata)) stop("Stratified Cox models are not currently supported")
    par <- coef(wide)
    inds <- check_inds(inds, length(par))
    inds0 <- check_inds0(inds0, inds, length(par))

    nmod <- nrow(inds)
    X <- check_X(X, length(par))
    t <- check_t(t)
    nval <- nrow(X)  # number of covariate values to evaluate the focus at
    if (isTRUE(sub=="auto"))
        sub <- fit_submodels(wide, inds)
    parsub <- get_parsub(sub, length(par), inds, inds0, gamma0, coef, wide)
    H0 <- basehaz(wide, centered=FALSE)
    fl <- get_focus_cox(focus, par=par, X=X, H0=H0, t=t)
    outn <- c("FIC", "rmse", "rmse.adj", "bias", "se")
    if (!is.null(parsub)) {
        parsub <- check_parsub(parsub, length(par), nmod)
        outn <- c(outn, "focus")
    }
    nout <- length(outn)  # number of outputs like FIC, rmse, rmse.adj,...
    res <- array(dim = c(length(t), nrow(X), nout, nmod))
    dimnames(res) <- list(tvals=t, vals=rownames(X), outn, mods=rownames(inds))
    for (i in 1:nmod){
        ficres <- fic_coxph_core(wide=wide, inds=inds[i,], inds0=inds0,
                                 gamma0=gamma0, focus=fl, X=X, t=t)
        focus_val <- if (!is.null(parsub)) fl$focus(par=parsub[i,], X=X, H0=attr(sub[[i]],"basehaz"), t=t) else NULL
        res[,,,i] <- abind(ficres, focus_val)
    }    
    if (tidy) {
        res <- tidy.array(res, dim2=3, ord=c("vals","tvals","mods"))
    }
    ## indices for wide and narrow models
    iwide <- apply(inds, 1, function(x)all(x==1))
    attr(res, "iwide") <- if (any(iwide)) which(iwide) else NULL
    inarr <- apply(inds, 1, function(x)all(x==inds0))
    attr(res, "inarr") <- if (any(inarr)) which(inarr) else NULL
    attr(res, "sub") <- sub
    attr(res, "parnames") <- get_parnames(par, inds)
    attr(res, "termnames") <- attr(inds, "termnames") 
    attr(res, "inds") <- inds

    class(res) <- c("fic",class(res))
    res
}

check_t <- function(t){
    if (is.null(t)) t <- 1 # e.g. for hazard ratio, which is indep of time
    if (!is.numeric(t)) stop("`t` must be numeric")
    if (any(t < 0)) stop("t must be a vector of survival times all > 0")
    t
}

fic_coxph_core <- function(wide,
                           inds,
                           inds0, 
                           gamma0 = 0,
                           focus = NULL, # list returned by get_focus_cox.
                           X=NULL,
                           t=NULL,
                           ...
                           ) 
{
    ## Notation in book 3.4
    ## Only handle built-in focuses for the moment 
    par <- coef(wide)
    n <- nrow(model.frame(wide))
    J <- solve(vcov(wide)) / n  # not divided by n here
    H0 <- basehaz(wide, centered=FALSE)

    if (is.null(t)) t <- 1 # dummy value 
    ncov <- nrow(X) # user-supplied number of covariate vectors to evaluate focus at 
    ntimes <- length(t)
    npar <- length(par)
    
    check_J(J, npar)
    inds <- check_indsinds0(inds,inds0)
    qq <- sum(inds0==0)  # maximum number of "extra" parameters
    gamma0 <- check_gamma0(gamma0, inds0)

    i0 <- which(inds0==1)
    dnhat <- sqrt(n)*(par[-i0] - gamma0) # deltahat / n
    J00 <- J[i0, i0, drop=FALSE]
    J10 <- J[-i0, i0, drop=FALSE]
    J01 <- J[i0,-i0, drop=FALSE]
    J11 <- J[-i0,-i0, drop=FALSE]
    invJ <- solve(J)
    Q <- invJ[-i0,-i0, drop=FALSE]   # using book notation.  called K in original code.

    focus_val <- focus
    focus_deriv <- focus$focus_deriv # npar x ntimes x ncov 
    dmudH0 <- focus$focus_dH # ntimes x ncov 

    survt <- model.frame(wide)[,attr(terms(model.frame(wide)), "response")]
    dead <- survt[,"status"]
    ti <- survt[,"time"] 
    deathtimes <- unique(ti[dead==1])
    H0u <- rbind(c(0,0), H0[H0$time %in% deathtimes,])
    pred <- linpred(wide)
    XZ <- model.matrix(wide)
    Yiu <- outer(ti, deathtimes, ">=")
    Gn0 <- colMeans(exp(pred) * Yiu)
    Gn1 <- t(XZ) %*% (exp(pred)*Yiu) / nrow(XZ)
    En <- t(t(Gn1) / Gn0) # npar x ndeaths
    tind <- outer(H0u$time[-1], t, "<=") # ndeaths x ntimes
    integrand <- t(En) * diff(H0u$hazard) # ndeaths x npar
    Fhat <- matrix(nrow=npar, ncol=ntimes)
    for (i in 1:npar){
        Fhat[i,] <- colSums(integrand[,i]*tind)
    }
    F0 <- Fhat[inds0==1,,drop=FALSE] # 
    F1 <- Fhat[inds0==0,,drop=FALSE]
    i0 <- which(inds0==1)
    i1 <- which(inds0==0)
    dmudbeta <- focus_deriv[i0,,,drop=FALSE]    # p x ntimes x ncov 
    dmudgamma <- focus_deriv[i1,,,drop=FALSE]  # q x ntimes x ncov
    pp <- length(i0)
    integrand <- H0u$hazard[-1] * diff(H0u$time) / Gn0
    intdHg <- colSums(integrand*tind)
    ## expand dimensions to match 
    dmudH0.rep <- array(dmudH0, dim=c(1, ntimes, ncov))[rep(1,pp),,,drop=FALSE]
    F0.rep <- array(F0, dim=c(pp, ntimes, ncov))
    dmF <- dmudbeta - dmudH0.rep*F0.rep
    dmFsq <- diag.array(tensor(tensor(dmF, J00, 1, 1), dmF, 3, 1))
    tau0sq <- dmudH0^2 * intdHg  +  dmFsq   # ntimes x ncov 
    ## J10 %*% solve(J00) is  q x p
    ## dmudbeta is p x ntimes x ncov
    ## F0 is p x ntimes,  F1 is q x ntimes 
    omega <- tensor(J10 %*% solve(J00), dmudbeta, 2, 1) - dmudgamma  # q x ntimes x ncov 
    JF <- J10 %*% solve(J00) %*% F0 - F1 # q x ntimes
    JF.rep <- array(JF, dim=c(qq, ntimes, ncov))
    ## dmudH0 is ntimes x ncov.  replicate to q x ntimes x ncov 
    dmudH0.rep <- array(dmudH0, dim=c(1, ntimes, ncov))[rep(1,qq),,,drop=FALSE]

    kappa <- JF.rep * dmudH0.rep # should be q x ntimes x ncov 
    ## dnhat is q vector, omega-kappa is q x ntimes x ncov
    ## this returns their "matrix product" with 1st dim of each summed over 
    psi.full <- tensor(omega - kappa, dnhat, 1, 1)
    
    ## similar to old FIC code but but replacing omega by omega - kappa 
    indsS <- inds[inds0==0]
    if (sum(indsS) > 0) { 
        Id <-  diag(rep(1,qq))
        Qinv <- solve(Q)  # q x q 
        pi.S <- matrix(Id[indsS*(seq_len(qq)),],ncol=qq) # qs x q
        Q.S <- solve(pi.S %*% Qinv %*% t(pi.S)) # qs x qs
        Q0.S <- t(pi.S) %*% Q.S %*% pi.S        # q x q
        G.S <- Q0.S %*% Qinv # called M.S in original code.  q x q
        ## G.S is q x q,  G.S %*% dnhat is q x 1
        psi.S <- tensor(omega - kappa, as.numeric(G.S %*% dnhat), 1, 1)

        ## unadjusted FIC estimates
        bias.S <- psi.full - psi.S # ntimes x ncov 
        ## array equivalent of
        ## FIC.S <- bias.S^2  +
        ##     2*diag(t(omega - kappa) %*% Q0.S %*% (omega - kappa))
        ## m x 1  +  diag of m x m
        ## where Q0.S is qxq, omega x kappa is q x ntimes x ncov  
        WQ0S <- tensor(omega - kappa, Q0.S, 1, 1)
        WQ0SW <- tensor(WQ0S, omega - kappa, 3, 1)
#        FIC.S <- bias.S^2  +  2*diag.array(WQ0SW)

        ## bias-adjusted estimates
        ## omega is q x ntimes x ncov,  IDI is q x q 
        IDI <- (Id - G.S) %*% (dnhat %*% t(dnhat) - Q) %*% t(Id - G.S) # q x q 
        sqbias2 <- tensor(tensor(omega - kappa, IDI, 1, 1), omega - kappa, 3, 1)
        sqbias2 <- array(diag.array(sqbias2), dim=c(ntimes, ncov))
        sqbias3 <- pmax(sqbias2, 0)
        bias.adj.S <- sign(bias.S) * sqrt(sqbias3)

        var.S <- tau0sq  +  diag.array(WQ0SW)
        mse.S <- sqbias2 + var.S
        WQW <- tensor(tensor(omega - kappa, Q, 1, 1), omega - kappa, 3, 1)
        FIC.S <- (mse.S - tau0sq + diag.array(WQW))
    } else { 
        ## Special case for null model with all extra parameters excluded
        bias.S <- bias.adj.S <- psi.full
        var.S <- tau0sq
        mse.S <- bias.S^2 + var.S
        FIC.S <- bias.S^2
    }    
    mse.adj.S <- var.S + bias.adj.S^2

    ## unadjusted mse
#    mse.S <- FIC.S + tau0sq - diag.array(WQW)
    res <- abind(
        FIC = FIC.S,
        rmse = sqrt_nowarning(mse.S/n),
        rmse.adj = sqrt(mse.adj.S/n),   # book p157
        bias = bias.adj.S/sqrt(n),
        se = sqrt(var.S/n),
        along=3
    )
    res
}

## Linear predictor from a Cox model
## note predict.coxph centres covariates around their means 

linpred <- function(mod, centered=FALSE){ 
    mm <- model.matrix(mod)
    if (centered) {
        xbar <- matrix(mod$means, nrow=nrow(mm), ncol=ncol(mm), byrow=TRUE)
        mm <- mm - xbar
    }
    as.numeric(mm %*% coef(mod))
}


## Generalization of matrix diagonal to 4D n x m x n x m arrays 
## Returns matrix with i,j element given by i,j,i,j element of array 
diag.array <- function(arr){
    dim1 <- dim(arr)[1]
    dim2 <- dim(arr)[2]
    eg <- expand.grid(seq(length=dim1), seq(length=dim2))
    inds <- as.matrix(cbind(eg, eg))
    array(arr[inds], dim=c(dim1,dim2))
}


get_focus_cox <- function(focus, focus_deriv=NULL, focus_dH=NULL, par=NULL, H0=NULL, X=NULL, t=NULL, ...){
    if (is.character(focus)) {
        if (!(focus %in% names(cox_focus_fns)))
            stop(sprintf("focus function `%s` not found", focus))
        fi <- match(focus, names(cox_focus_fns))
        args_extra <- list(...)
        focus <- cox_focus_fns[[fi]]$focus
        deriv_args <- c(par=list(par), H0=list(H0), X=list(X), t=list(t), args_extra)
        focus_deriv <- do.call(cox_focus_fns[[fi]]$deriv, deriv_args)
        focus_dH <- do.call(cox_focus_fns[[fi]]$dH, deriv_args)
        if (!("X" %in% names(formals(focus))))
            formals(focus) <- c(formals(focus), alist(X=))
        if (!("t" %in% names(formals(focus))))
            formals(focus) <- c(formals(focus), alist(t=))
        res <- list(focus=focus, focus_deriv=focus_deriv, focus_dH=focus_dH)
    } else if (is.list(focus)){
        args <- list(par, H0, X, t)
        res <- list(
            focus = focus$focus,
            focus_deriv = focus$focus_deriv(par, H0, X, t),
            focus_dH = focus$focus_dH(par, H0, X, t))
    }
    res
    ## add unused X and t argument to function if it doesn't already have one
}

## HRs between given X and 0
## by definition, these don't depend on time
## ntimes = 1 usually

cox_hr <- function(par, H0, X, t) {
    ntimes <- length(t)
    ncov <- nrow(X)
    matrix(exp(X %*% par), nrow=ntimes, ncol=ncov, byrow=TRUE)
}

cox_hr_deriv <- function(par, H0, X, t){
    npar <- length(par)
    ntimes <- length(t) 
    ncov <- nrow(X)
    e <- as.vector(exp(X %*% par))
    array(t(X * e), dim=c(npar, ntimes, ncov))
}

cox_hr_dH <- function(par, H0, X, t){
    ntimes <- length(t)
    ncov <- nrow(X)
    array(0, dim=c(ntimes, ncov))
}

##' Interpolate cumulative hazard function from a fitted Cox model
##'
##' Returns the baseline cumulative hazard, at the requested times,
##' from a Cox model fitted by \code{\link[survival]{coxph}}.  Linear
##' interpolation is used, assuming the hazard is piecewise constant,
##' thus the cumulative hazard is piecewise linear.
##'
##' @param H0 output from \code{\link[survival]{basehaz}}, containing estimates of the baseline cumulative hazard at a series of times.
##'
##' @param t vector of times for which cumulative hazard estimates are required.
##'
##' @return Fitted cumulative hazard at \code{t}.
##'
##' @details This does not extrapolate.  If \code{t} is outside the observed event times, then \code{NA} will be returned.
##'
##' @examples
##'
##' library(survival)
##' wide <- coxph(Surv(years, death==1) ~ sex + thick_centred +
##'               infilt + epith + ulcer + depth + age, data=melanoma)
##' basehaz(wide)
##' get_H0(basehaz(wide), c(0,1,5,10,100))
##' 
##' 
##' @export
##' 
get_H0 <- function(H0, t){
    H0 <- rbind(c(0, 0), H0) # cum haz always 0 at time 0 
    approx(H0$time, H0$hazard, xout=t)$y
}

cox_cumhaz <- function(par, H0, X, t) {
    H0t <- get_H0(H0, t)
    outer(H0t, exp(X %*% par)[,1])  # ntimes x ncov 
}

cox_cumhaz_deriv <- function(par, H0, X, t) {
    H0t <- get_H0(H0, t)
    ncov <- nrow(X)
    npar <- length(par)
    ntimes <- length(t) 
    ## X is ncov x npar
    Xrep <- array(t(X), dim=c(npar, 1, ncov))[,rep(1,ntimes),,drop=FALSE] # npar x ntimes x ncov 
    out <- outer(H0t, exp(X %*% par)) # ntimes x ncov x 1
    out.rep <- array(out, dim=c(1,ntimes,ncov))[rep(1,npar),,,drop=FALSE] # npar x ntimes x ncov
    Xrep * out.rep
}

cox_cumhaz_dH <- function(par, H0, X, t) {
    ntimes <- length(t)
    ncov <- nrow(X)
    array(1, dim=c(ntimes, ncov))
}

cox_survival <- function(par, H0, X, t) {
    exp(-cox_cumhaz(par, H0, X, t))   # ntimes x ncov 
}

cox_survival_deriv <- function(par, H0, X, t) {
    ncov <- nrow(X)
    npar <- length(par)
    ntimes <- length(t)
    surv <- cox_survival(par, H0, X, t) # ntimes x ncov 
    surv.rep <- array(surv, dim=c(1, ntimes, ncov))[rep(1,npar),,,drop=FALSE]
    - surv.rep * cox_cumhaz_deriv(par, H0, X, t) # both terms:  npar x ntimes x ncov 
}

cox_survival_dH <- function(par, H0, X, t) {
    ntimes <- length(t)
    hr <- t(exp(X %*% par))  # ncov x 1 transposed to 1 x ncov
    hr.rep <- hr[rep(1,ntimes),,drop=FALSE]  # ntimes x ncov
    - cox_survival(par, H0, X, t) * hr.rep      # ntimes x ncov
}

cox_focus_fns <- list(
    "hr" = list(focus        = cox_hr,
                deriv  = cox_hr_deriv,
                dH     = cox_hr_dH),
    "cumhaz" = list(focus    = cox_cumhaz,
                deriv  = cox_cumhaz_deriv,
                dH     = cox_cumhaz_dH),
    "survival" = list(focus   = cox_survival,
                deriv  = cox_survival_deriv,
                dH     = cox_survival_dH)
)
