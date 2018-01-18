##' Focused Information Criterion: generic function
##' 
##' 
##' @rdname fic
##' @export fic
fic <- function(x,...) UseMethod("fic")
  


##' Focused Information Criterion: default method
##'
##' Focused information criterion for general models.  These methods estimate the resulting bias and variance of estimates of a quantity of interest (the "focus") when parameters are excluded from a "wide" model that is assumed to encompass all plausible models.
##'
##' @aliases FIC
##' 
##' @param par Vector of maximum likelihood estimates from the wide model
##' 
##' @param J Information matrix from the wide model, evaluated at the maximum likelihood estimates and divided by \code{n}.
##'
##' @param inds Vector of 0s and 1s of same length as \code{par}, with 1s in the positions where the parameters of the wide model are included in the submodel to be assessed, and 0s elsewhere.
##'
##' @param inds0 Specification of narrow model, in the same format as \code{inds}.  TODO error checking
##' 
##' @param n Number of observations in the data used to fit the wide model.  TODO is this really needed, or does it cancel?
##'
##' @param focus An R function whose first argument is a vector of parameters; and which returns one scalar: the focus quantity of interest.  Not required if \code{focus_deriv} is specified.
##'
##' Alternatively, this can be a character string naming a built-in focus function supplied by the \pkg{fic} package.  See \code{\link{focus_fns}}. 
##'
##' @param focus_deriv Vector of partial derivatives of the focus function with respect to the parameters in the wide model.  If not supplied, this is computed from the function supplied in \code{focus}, using numerical differentiation.
##'
##' @param parsub Vector of maximum likelihood estimates from the submodel.  
##' Only required to return the estimate of the focus quantity alongside 
##' the model assessment statistics for the submodel. If omitted, the estimate is omitted.
##'
##' @param gamma0 Vector of special values taken by gamma in the narrow model.  
##' Defaults to all 0, as in covariate selection where coefficients are fixed to 0 
##' Either a scalar, assumed to be the same for all elements of gamma, or a vector or the same length as gamma.
##' Ignored unless \code{parsub} is specified.
##'
##' @param \dots Other arguments to the focus function can be supplied here.
##'
##' These might include the covariate values at which the focus is to be evaluated.  For the built-in focus functions, this argument is named \code{X} and can either be a vector of length \eqn{p}, or a \eqn{n x p} matrix, where \eqn{p} is the number of covariate effects, and \eqn{n} is the number of alternative sets of covariate values at which the focus function, and hence FIC, is to be evaluated.  (TODO not implemented vectorisation yet)
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
fic.default <- function(par, # estimates in wide model 
                J, #  info matrix in wide model divided by n
                inds, #  indicator of which terms to include in submodel (was called "variables")
                inds0, 
                n,
                focus=NULL,
                focus_deriv=NULL,
                parsub=NULL, # estimates in submodel
                gamma0=0,
                ...
                ) 
{
    pp <- sum(inds0)
    qq <- sum(inds0==0)  # maximum number of "extra" covariates
    i0 <- which(inds0==1)
    indsS <- inds[inds0==0]
    if (any(inds[inds0==1] != 1)){
      dodgy_inds <- paste(which(inds[inds0==1] != 1), collapse=",")
      warning("Submodel excludes parameters in the narrow model, in position ", dodgy_inds, ". Carrying on and including them. ")
      inds[inds0==1] <- 1
    }
    
    deltahat <- sqrt(n)*par[-i0]
    J00 <- J[i0, i0]
    J10 <- J[-i0, i0]
    J01 <- J[i0,-i0]
    J11 <- J[-i0,-i0]
    invJ <- solve(J)
    Q <- invJ[-i0,-i0]   # using book notation.  called K in original code.
   
    # Handle built-in focus functions
    fl <- get_focus(focus, focus_deriv, par, ...)
    focus <- fl$focus
    focus_deriv <- fl$focus_deriv
    
    if(is.null(focus_deriv))
        focus_deriv <- numDeriv::grad(func=focus, x=par)
    dmudtheta <- focus_deriv[i0]
    tau0sq <- t(dmudtheta) %*% solve(J00) %*% dmudtheta
    dmudgamma <- focus_deriv[pp + seq_len(qq)]
    omega <- J10 %*% solve(J00) %*% dmudtheta - dmudgamma
    psi.full <- t(omega) %*% deltahat

    if (sum(indsS) > 0) { 
        Id <-  diag(rep(1,qq))
        Qinv <- solve(Q)
        pi.S <- matrix(Id[indsS*(seq_len(qq)),],ncol=qq)

        Q.S <- solve(pi.S %*% Qinv %*% t(pi.S))
        Q0.S <- t(pi.S) %*% Q.S %*% pi.S
        G.S <- Q0.S %*% Qinv # called M.S in original code
        psi.S <- t(omega)%*% G.S %*% deltahat

        omega.S <- pi.S %*% omega

        ## basic FIC estimates (section 6.1) 
        bias.S <- psi.full - psi.S
        FIC.S <- bias.S^2  +  2*t(omega.S) %*% Q.S %*% omega.S
        
        ## bias-adjusted estimates (section 6.4)
        sqbias2 <- t(omega) %*% (Id - G.S) %*% (deltahat %*% t(deltahat) - Q) %*% (Id - G.S) %*% omega # \hat{sqb2}(S) on p152 
        sqbias3 <- max(sqbias2, 0) # \hat{sqb2}(S) on p152
        bias.adj.S <- sign(bias.S) * sqrt(sqbias3)
        ## using Q0.S here as in book, rather than G.S (called M.S in code) as in original code.  is this right?
        var.S <- tau0sq  +  t(omega) %*% Q0.S %*% omega 
    } else { 
        ## Special case for null model with all extra parameters excluded
        bias.S <- bias.adj.S <- psi.full
        var.S <- tau0sq
        FIC.S <- bias.S^2
    }    
    mse.adj.S <- var.S + bias.adj.S^2

    ## unadjusted mse
    mse.S <- FIC.S + tau0sq - t(omega) %*% Q %*% omega # book p150, eq 6.7 
    
    res <- c(FIC = FIC.S,
             rmse = sqrt(mse.S / n),
             rmse.adj = sqrt(mse.adj.S / n),   # book p157
             bias = bias.S / sqrt(n),
             bias.adj = bias.adj.S / sqrt(n),
             se = sqrt(var.S / n)
             )
    if (!is.null(parsub)) 
      res["focus"] <- focus_sub(focus=focus, parsub=parsub, 
                                inds=inds, inds0=inds0, gamma0=gamma0, ...)
    res
}

FIC <- fic


##' Calculate the focus estimate under a submodel of the wide model
##' 
##' @inheritParams fic
##' 
focus_sub <- function(focus, parsub, inds, inds0, gamma0=NULL, ...){
    pp <- sum(inds0==1)
    qq <- sum(inds0==0) # maximum number of "extra" covariates
    if (is.null(gamma0)) gamma0 <- 0
    gamma0 <- rep(gamma0, length.out = qq)
    ests_long <- numeric(pp + qq)
    ests_long[inds==1] <- parsub
    ## indices of gamma in full parameters
    gi <- which(inds0==0)
    ## indices of gamma not in submodel
    ngi <- which(inds0==0 & inds==0)
    ests_long[inds0==0][gi %in% ngi] <- gamma0[gi %in% ngi]
    focus(ests_long, ...)
}
