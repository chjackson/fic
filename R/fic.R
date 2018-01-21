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
##' @param X second argument of the focus function.  This will typically be a vector of covariate values.  
##' 
##' This can be a vector with \code{nc} elements, e.g. one for each covariate coefficient in the wide model.
##' 
##' Alternatively this can be a matrix with \code{nf} rows and \code{nc} columns. 
##' The FIC calculation is then "vectorised".  With one call of the \code{fic} function, the   
##' FIC and related statistics are calculated \code{nf} times, with 
##' the focus function evaluated at each of \code{nf} alternative 
##' covariate values. 
##' 
##' Not yet tested in situations other than supplying covariate values. 
##' 
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
fic.default <- function(par, # estimates in wide model 
                J, #  info matrix in wide model divided by n
                inds, #  indicator of which terms to include in submodel (was called "variables")
                inds0, 
                n,
                focus=NULL,
                focus_deriv=NULL,
                parsub=NULL, # estimates in submodel
                gamma0=0,
                X=NULL, # second argument of focus function, e.g. covariate values 
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
    
    ## TODO check dims of X 
    ## Should be a vector of length ncoefs, or nvals x ncoefs matrix
    if(is.null(focus_deriv)){
      if(is.vector(X)) X <- matrix(X, nrow=1)
      ncols <- if(is.null(X)) 1 else nrow(X)
      focus_deriv <- matrix(nrow=length(par), ncol=ncols)
      for (i in 1:ncols){
        if (is.null(X))
          focus_deriv[,i] <- numDeriv::grad(func=focus, x=par) 
        else
          focus_deriv[,i] <- numDeriv::grad(func=focus, x=par, X=X[i,]) 
      }
    }
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
    colnames(res) <- c("FIC.S", "rmse", "rmse.adj", "bias", "bias.adj", "se")

    if (!is.null(parsub)) {
      fval <- focus_sub(focus=focus, parsub=parsub, inds=inds, inds0=inds0, X=X, gamma0=gamma0, ...)
      res <- cbind(res, focus=as.numeric(fval))
    }
    res
}

FIC <- fic


##' Calculate the focus estimate under a submodel of the wide model
##' 
##' @inheritParams fic
##' 
focus_sub <- function(focus, parsub, inds, inds0, X=NULL, gamma0=NULL, ...){
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
    if ("X" %in% names(formals(focus)))
      res <- focus(ests_long, X=X, ...)
    else
      res <- focus(ests_long, ...)
    res
}
