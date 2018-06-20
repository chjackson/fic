##' @aliases fic_core
##' 
##' @inheritParams fic_multi
##' 
##' @rdname fic_multi
##' @export
fic_core <- function(
                       par,
                       J, 
                       inds,
                       inds0, 
                       gamma0 = 0,
                       n,
                       focus = NULL,
                       focus_deriv = NULL,
                       X = NULL, 
                       ...
                       ) 
{
    npar <- length(par)
    if (!is.numeric(inds0)) 
        stop("`inds0` must be numeric")
    if (!is.vector(inds0))
        stop("`inds0` must be a vector")
    if (length(inds0) != npar)
        stop(sprintf("`inds0` of length %d, but model has %d parameters.\nLength of `inds0` must match number of parameters", length(inds0), npar))
    if (!is.matrix(J) || (nrow(J)!=ncol(J)))
        stop("`J` must be a square matrix")
    if (!is.numeric(J))
        stop("`J` must be numeric")
    if (nrow(J) != npar)
        stop(sprintf("`J` has %d rows and columns, but there are %d parameters.\n`J` must have number of rows and columns equal to the number of parameters", nrow(J), npar))

    pp <- sum(inds0)
    qq <- sum(inds0==0)  # maximum number of "extra" covariates
    
    if ((length(gamma0) != 1) && (length(gamma0) != qq))
        stop(sprintf("`gamma0` of length %d, but `inds0` has %d entries which are zero.\nLength of gamma0 must either be 1 or match the number of entries of `inds0` which are zero", length(gamma0), qq))

    i0 <- which(inds0==1)
    indsS <- inds[inds0==0]
    if (any(inds[inds0==1] != 1)){
      dodgy_inds <- paste(which(inds[inds0==1] != 1), collapse=",")
      warning("Submodel excludes parameters in the narrow model, in position ", dodgy_inds, ". Carrying on and including them. ")
      inds[inds0==1] <- 1
    }

    gamma0 <- rep(gamma0, length.out = qq) # TODO error checking 
    deltahat <- sqrt(n)*(par[-i0] - gamma0)
    J00 <- J[i0, i0]
    J10 <- J[-i0, i0]
    J01 <- J[i0,-i0]
    J11 <- J[-i0,-i0]
    invJ <- solve(J)
    Q <- invJ[-i0,-i0]   # using book notation.  called K in original code.
    
    dmudtheta <- focus_deriv[i0,]
    tau0sq <- diag(t(dmudtheta) %*% solve(J00) %*% dmudtheta)
    dmudgamma <- focus_deriv[pp + seq_len(qq),]

    omega <- J10 %*% solve(J00) %*% dmudtheta - dmudgamma # q x m 
    psi.full <- t(omega) %*% deltahat # m x 1

    if (sum(indsS) > 0) { 
        Id <-  diag(rep(1,qq))
        Qinv <- solve(Q)  # q x q 
        pi.S <- matrix(Id[indsS*(seq_len(qq)),],ncol=qq) # qs x q

        Q.S <- solve(pi.S %*% Qinv %*% t(pi.S)) # qs x qs
        Q0.S <- t(pi.S) %*% Q.S %*% pi.S        # q x q
        G.S <- Q0.S %*% Qinv # called M.S in original code.  q x q
        psi.S <- t(omega)%*% G.S %*% deltahat  # m x 1 

        omega.S <- pi.S %*% omega              # qs x m

        ## basic FIC estimates (section 6.1) 
        bias.S <- psi.full - psi.S # m x 1 
        FIC.S <- bias.S^2  +  2*diag(t(omega.S) %*% Q.S %*% omega.S)  # m x 1  +  diag of m x m 
        
        ## todo test all diagonalisations

        ## bias-adjusted estimates (section 6.4)
        sqbias2 <- t(omega) %*% (Id - G.S) %*% (deltahat %*% t(deltahat) - Q) %*% (Id - G.S) %*% omega # \hat{sqb2}(S) on p152 
        sqbias2 <- diag(sqbias2)
        sqbias3 <- pmax(sqbias2, 0) # \hat{sqb2}(S) on p152
        bias.adj.S <- sign(bias.S) * sqrt(sqbias3)
        ## using Q0.S here as in book, rather than G.S (called M.S in code) as in original code.  is this right?
        var.S <- tau0sq  +  diag(t(omega) %*% Q0.S %*% omega)
    } else { 
        ## Special case for null model with all extra parameters excluded
        bias.S <- bias.adj.S <- psi.full
        var.S <- tau0sq
        FIC.S <- bias.S^2
    }    
    mse.adj.S <- var.S + bias.adj.S^2

    ## unadjusted mse
    mse.S <- FIC.S + tau0sq - diag(t(omega) %*% Q %*% omega) # book p150, eq 6.7 
    res <- cbind(
             FIC      = FIC.S,
             rmse     = sqrt(mse.S / n),
             rmse.adj = sqrt(mse.adj.S / n),   # book p157
             bias     = bias.S / sqrt(n),
             bias.adj = bias.adj.S / sqrt(n),
             se       = sqrt(var.S / n)
             )
    colnames(res) <- c("FIC", "rmse", "rmse.adj", "bias", "bias.adj", "se")
        
    res
}


check_inds <- function(inds, npar){
    if (!is.numeric(inds)) stop("`inds` must be numeric")
    if (!(is.vector(inds) || is.matrix(inds)))
        stop("`inds` must be a vector or a matrix")
    if (is.vector(inds)){
        if (length(inds) != npar)
            stop(sprintf("`inds` of length %d, but model has %d parameters.\nLength of vector `inds` must match number of parameters", length(inds), npar))
        inds <- matrix(inds, nrow=1)
    } else {
        if (ncol(inds) != npar)
            stop(sprintf("`inds` has %d columns, but model has %d parameters.\nNumber of columns of matrix `inds` must match number of parameters", ncol(inds), npar))
    }
    if (is.null(rownames(inds))) rownames(inds) <- seq_len(nrow(inds))
    inds
}


##' Focused Information Criterion: core calculation functions
##'
##' Core FIC calculation functions underlying the user interface in \code{\link{fic}}.
##' \code{fic_core} just handles one submodel, while \code{fic_multi} can assess multiple submodels of the same wide model.  For \code{fic_multi}, \code{\link{inds}} and \code{\link{parsub}} can be matrices with one row per submodel, while for \code{fic_core} they must be vectors.
##' 
##' @param par Vector of maximum likelihood estimates from the wide model
##' 
##' @param J Information matrix from the wide model, evaluated at the maximum likelihood estimates and divided by \code{n}.
##'
##' @param inds Vector of 0s and 1s of same length as \code{par}, with 1s in the positions where the parameters of the wide model are included in the submodel to be assessed, and 0s elsewhere.
##'
##' @param n Number of observations in the data used to fit the wide model.  TODO is this really needed for the MSE, or does it cancel?  Only of academic interest to compute FIC/error for the sqrt n transformed quantity? 
##'
##' @param parsub Vector of maximum likelihood estimates from the submodel.  
##' Only required to return the estimate of the focus quantity alongside 
##' the model assessment statistics for the submodel. If omitted, the estimate is omitted.
##'@param focus_deriv Vector of partial derivatives of the focus function with respect to the parameters in the wide model.  This is required by \code{fic.core}.
##' 
##' If there are multiple submodels, this should be a matrix with number of rows equal to the number of submodels, and number of columns equal to the number of parameters in the wide model.  If there is a single submodel, this should be a vector with number of columns equal to the number of parameters in the wide model.
##'
##' This should take the value given by \code{gamma0} in positions where the parameter is excluded from the submodel.
##' For example, coefficients of covariates should take the value 0 if the covariate is excluded. 
##'
##' @inheritParams fic.default
##'
##' @return See \code{\link{fic}} 
##' 
##' @rdname fic_multi
##'
##' @export
fic_multi <- function(
                       par,
                       J, 
                       inds,
                       inds0, 
                       gamma0 = 0,
                       n,
                       focus = NULL,
                       focus_deriv = NULL,
                       X = NULL, 
                       parsub = NULL,
                       ...
                      )

{
    npar <- length(par)
    inds <- check_inds(inds, npar)
    nmod <- nrow(inds)

    ## TODO error check for dims of X (ncoef <= npar?) 
    if (!is.null(X)){
        if (!is.numeric(X)) stop("`X` must be numeric")
        if (!(is.vector(X) || is.matrix(X))) stop("X must be a vector or a matrix")
        if (is.vector(X)) {
            X <- matrix(X, nrow=1)
            nval <- 1
            ncoef <- length(X)
        } else {
            nval <- nrow(X)
            ncoef <- ncol(X)
        }
        if (is.null(rownames(X))) rownames(X) <- seq_len(nrow(X))
    } else nval <- 1

    outn <- c("FIC", "rmse", "rmse.adj", "bias", "bias.adj", "se")

    if (!is.null(parsub)) {
        ## TODO more error checks
        if (!is.numeric(parsub)) stop("`parsub` must be numeric")
        outn <- c(outn, "focus")
    }
   
    # Handle built-in focus functions
    fl <- get_focus(focus, focus_deriv, par, X, ...)
    focus <- fl$focus
    focus_deriv <- fl$focus_deriv
    if (!is.function(focus)) stop("`focus` must be a function")
    if(is.null(focus_deriv)){
      ncols <- if(is.null(X)) 1 else nrow(X)
      focus_deriv <- matrix(nrow=npar, ncol=ncols)
      for (i in 1:ncols){
        if (is.null(X))
          focus_deriv[,i] <- numDeriv::grad(func=focus, x=par) 
        else
          focus_deriv[,i] <- numDeriv::grad(func=focus, x=par, X=X[i,]) 
      }
    }

    nout <- length(outn)
    res <- array(dim = c(nval, nout, nmod))
    dimnames(res) <- list(vals=rownames(X), outn, mods=rownames(inds))
    for (i in 1:nmod){
        ficres <- fic_core(par=par, J=J, inds=inds[i,], inds0=inds0,
                             gamma0=gamma0, n=n, focus=focus, focus_deriv=focus_deriv,
                             X=X, ...)
        
        fres <- if (!is.null(parsub)) focus(par=parsub[i,], X=X, ...) else NULL
        res[,,i] <- cbind(ficres, fres)
    }
    res
}

get_fns <- function(fns){
    args <- c("coef","nobs","vcov","AIC","BIC")
    default_fns <- list(
        coef = coef,
        nobs = nobs,
        vcov = vcov,
        AIC = AIC,
        BIC = BIC
    )
    ret <- default_fns
### TODO error checking. must be list of functions, names must be in args 
    if (!is.null(fns)){
        for (i in names(fns))
            ret[i] <- fns[i]
    }
    ret
}

get_parsub <- function(sub, npar, inds, inds0, gamma0, coef_fn){
    if (is.null(sub))
        parsub <- NULL
    else {
        ## todo error checking for sub of correct form
        nmod <- length(sub)
        parsub <- array(dim=c(nmod, npar))
        for (i in 1:nmod){
            parsub[i,inds0==0] <- gamma0
            parsub[i,inds[i,]==1] <- coef_fn(sub[[i]])
        }
    }
    parsub
}

# Not sure the package really needs to return these
get_ics <- function(sub, fns){
    if (!is.list(sub))
        sub <- list(sub)
    nmod <- length(sub)
    aics <- sapply(sub, fns$AIC)
    bics <- sapply(sub, fns$BIC)
}

##' Focused Information Criterion: main user interface
##'
##' Focused information criterion for general models.  These methods estimate the bias and variance of estimates of a quantity of interest (the "focus") when parameters are excluded from an "wide" model that is assumed to generate the data.
##'
##' @aliases FIC
##'
##' @param wide Fitted model object containing the wide model.
##'
##' @param inds matrix or vector of indicators for which parameters are included in the submodel or submodels to be assessed.
##'
##' A matrix should be supplied if multiple submodels are to be assessed.  This should have number of rows equal to the number of submodels to be assessed, and number of columns equal to the total number of parameters in the wide model.  It contains 1s in the positions where the parameter is included in the submodel, and 0s in positions where the parameter is excluded.  This should always be 1 in the positions defining the narrow model, as specified in `inds0`.
##'
##' @param inds0 Specification of narrow model, in the same format as \code{inds}.  TODO error checking
##'
##' @param gamma0 Vector of special values taken by the parameters \eqn{gamma} which define the narrow model.
##' 
##' This defaults to 0, as in covariate selection, where excluded coefficients are fixed to 0. 
##'
##' This should either be a scalar, assumed to be the same for all parameters fixed in the narrow model, or a vector of length equal to the number of parameters from the wide model which are fixed in the narrow model, that is, the number of entries of inds0 which are zero.
##' 
##'
##' @param focus An R function with:
##'
##' \itemize{
##' \item first argument named \code{par}, denoting a vector of parameters
##'
##' \item an optional second argument named \code{X}, typically denoting covariate values.  The required format is documented below. 
##' }
##' 
##' The function should return the focus quantity of interest.  If \code{X} is supplied and has multiple rows, then \code{focus} should return a vector giving the focus for each row of \code{X}.  Otherwise \code{focus} should return a scalar giving the focus value at \code{par}.
##'
##' Not required if \code{focus_deriv} is specified.
##'
##' Alternatively, \code{focus} can be a character string naming a built-in focus function supplied by the \pkg{fic} package.  See \code{\link{focus_fns}}. 
##'
##' @param focus_deriv Vector of partial derivatives of the focus function with respect to the parameters in the wide model.  This is not usually needed, as it can generally be computed automatically and accurately from the function supplied in \code{focus}, using numerical differentiation.
##'
##' @param X Second argument of the focus function, typically giving covariate values defining the focus. This can either be a matrix or a vector.
##'
##' This is optional.  If not supplied, \code{focus} should return a scalar giving the focus value at parameters \code{par}.
##'
##' To compute focused model comparison statistics for the same focus function evaluated at multiple covariate values, \code{X} should be a matrix, with number of columns equal to the number of parameters in the wide model, and number of rows equal to the number of alternative covariate values.
##'
##' For a typical regression model, the first parameter will denote an intercept, so the first value of \code{X} should be 1, and the remaining values should correspond to covariates whose coefficients form parameters of the wide model.  See the examples in the vignette.
##'
##' If just one covariate value is needed, then \code{X} can be a vector of length equal to the number of parameters in the wide model. 
##'
##' @param sub List of fitted model objects for each submodel to be assessed.  This is optional, and only required if you want the estimate of the focus function under each submodel to be included in the results.
##'
##' @param fns Named list of functions to extract the quantities from the fitted model object that are required for the FIC calculation.  By default this is
##'
##' \code{list(coef=coef, nobs=nobs, vcov=vcov)}
##'
##' Suppose the fitted model object is called \code{mod}.  This default list assumes that
##'
##' \itemize{
##' \item \code{coef(mod)} returns the vector of parameter estimates,
##'
##' \item \code{nobs(mod)} returns the number of observations used in the model fit,
##'
##' \item \code{vcov(mod)} returns the covariance matrix for the parameter estimates.
##' }
##'
##' If one or more of these assumptions is not true, then the defaults can be changed.
##' For example, suppose the functions \code{coef()}, \code{nobs()} and \code{vcov()} are
##' not understood (or return something different) for your class of model objects, but
##' the parameters are stored in \code{mod$estimates},
##' the number of observations is in \code{mod$data$nobs}, and
##' the covariance matrix is in \code{mod$vcov}, 
##' 
##' then the \code{fns} argument should be set to 
##'
##' \code{list(
##'      coef = function(x){x$estimates},
##'      nobs = function(x){x$data$nobs},
##'      vcov = function(x){x$cov}
##'      )}
##'
##' If less than three components are specified in \code{fns}, then the missing components are assumed to take their default values.
##'
##' @param tidy If \code{TRUE} return the results as a data frame.  If \code{FALSE}, return the results as a three-dimensional array, with dimensions indexed by the submodels, result statistics and \code{X} values respectively.
##' 
##' @param \dots Other arguments to the focus function can be supplied here.  Currently no examples of this. 
##'
##' @return A vector containing the following components, describing characteristics of the defined submodel (references in Chapter 6 of Claeskens and Hjort, 2008)
##'
##' \item{FIC}{The focused information criterion (equation 6.1). }
##'
##' \item{rmse}{The root mean square error of the estimate of the focus quantity.  Defined on page 150 (equation 6.7).}
##'
##' \item{rmse.adj}{The root mean square error, with an adjustment to avoid negative squared bias.  Defined on page 157 as the sum of the variance and the squared adjusted bias.}
##'
##' \item{bias}{The estimated bias of the focus quantity (unadjusted).  Defined on page 157.}
##'
##' \item{bias.adj}{The estimated bias of the focus quantity (adjusted to avoid negative squared bias).  This is defined as the square root of the quantity "sqb3(S)", page 152, multiplied by the sign of the unadjusted bias, and divided by the square root of the sample size. }
##'
##' \item{se}{The estimated standard error (root variance) of the focus quantity.  Defined on page 157.}
##'
##' 
##' @references Claeskens, G., & Hjort, N. L. (2008). Model selection and model averaging (Vol. 330). Cambridge: Cambridge University Press.
##'
##' Claeskens, G., & Hjort, N. L. (2003). The focused information criterion. Journal of the American Statistical Association, 98(464), 900-916.
##'
##' @keywords models
##' 
##' @examples TODO
##'
##' @rdname fic
##' @export
fic.default <- function(wide, inds, inds0, gamma0=0, 
                      focus=NULL, focus_deriv=NULL,
                      X=NULL, sub=NULL, fns=NULL,
                      tidy=TRUE, ...){
    fns <- get_fns(fns)
    par <- fns$coef(wide)
    n <- fns$nobs(wide)
    J <- solve(fns$vcov(wide)) / n
    inds <- check_inds(inds, length(par))
    parsub <- get_parsub(sub, length(par), inds, inds0, gamma0, fns$coef)
    res <- fic_multi(par=par, J=J, inds=inds, inds0=inds0, gamma0=gamma0, n=n, 
              focus=focus, focus_deriv=focus_deriv, 
              parsub=parsub, X=X, ...)
    if (tidy){
        res1 <- apply(res, 2, as.data.frame.table)
        res <- do.call("data.frame", lapply(res1, function(x)x$Freq))
        res <- cbind(res1[[1]][,1:2], res)
    }
    res
}

##' @rdname fic
##' @export fic
fic <- function(x,...) UseMethod("fic")

FIC <- fic
fic_multi <- fic_multi
fic_core <- fic_core
