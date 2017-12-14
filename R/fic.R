##' Focused Information Criterion: general function
##'
##' Focused information criterion for general models.  These methods estimate the resulting bias and variance of estimates of a quantity of interest (the "focus") when parameters are excluded from a "wide" model that is assumed to encompass all plausible models.
##'
##' @aliases FIC
##' 
##' @param ests Vector of maximum likelihood estimates from the wide model
##'
##' @param J Information matrix from the wide model, evaluated at the maximum likelihood estimates and divided by \code{n}.
##'
##' @param inds Vector of 0s and 1s of length \code{length(ests) - pp}, with 1s in the positions where the parameters of the wide model are included in the submodel, and 0s in the positions where the parameters of the wide model are excluded from the submodel.
##'
##' @param pp Number of parameters which we would always include in any submodel.  The corresponding parameters are assumed to be the first \code{pp} parameters included in \code{ests}.
##'
##' @param n Number of observations in the data used to fit the wide model.
##'
##' @param focus An R function whose first argument is a vector of parameters; and which returns one scalar: the focus quantity of interest.  Not required if \code{focus_deriv} is specified.
##'
##' Alternatively, this can be a character string naming a built-in focus function supplied by the \pkg{fic} package.  See \code{\link{focus_fns}}. 
##'
##' @param focus_deriv Vector of partial derivatives of the focus function with respect to the parameters in the wide model.  If not supplied, this is computed from the function supplied in \code{focus}, using numerical differentiation.
##'
##' @param \dots Other arguments to the focus function can be supplied here.
##'
##' These might include the covariate values at which the focus is to be evaluated.  For the built-in focus functions, this argument is named \code{X} and can either be a vector of length \eqn{p}, or a \eqn{n x p} matrix, where \eqn{p} is the number of covariate effects, and \eqn{n} is the number of alternative sets of covariate values at which the focus function, and hence FIC, is to be evaluated.  (TODO not implemented vectorisation yet)
##'
##' @return A vector containing the following components, describing characteristics of the defined submodel (references in Chapter 6 of Claeskens and Hjort, 2008)
##'
##' \item{FIC}{The focused information criterion (equation 6.1). }
##'
##' \item{rmse}{The root mean square error of the estimate of the focus quantity.  Defined on page 157 as the sum of the variance and the squared adjusted bias.}
##'
##' \item{bias}{The estimated bias of the focus quantity (unadjusted).  Defined on page 157.}
##'
##' \item{bias.adj}{The estimated bias of the focus quantity (adjusted to avoid negative squared bias).  This is defined as the square root of the quantity "sqb3(S)", page 152, multiplied by the sign of the unadjusted bias, and divided by the square root of the sample size. }
##'
##' \item{se.adj}{The estimated standard error (root variance) of the focus quantity (adjusted).  Defined on page 157.}
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
##' @export
fic <- function(ests, # estimates in wide model 
                J, #  info matrix in wide model divided by n
                inds, #  indicator of which terms to include in submodel (was called "variables")
                pp, 
                n,
                focus=NULL,
                focus_deriv=NULL,
                ...
                ) 
{
    deltahat <- sqrt(n)*ests[-(1:pp)]
    J00 <- J[1:pp, 1:pp]
    J10 <- J[-(1:pp), 1:pp]
    J01 <- J[1:pp,-(1:pp)]
    J11 <- J[-(1:pp),-(1:pp)]
    invJ <- solve(J)
    Q <- invJ[-(1:pp),-(1:pp)]   # using book notation.  called K in original code.
    qq <- length(inds) # maximum number of "extra" covariates

    # Handle built-in focus functions
    fl <- get_focus(focus, focus_deriv, ests, ...)
    focus <- fl$focus
    focus_deriv <- fl$focus_deriv
    
    if(is.null(focus_deriv))
        focus_deriv <- numDeriv::grad(func=focus, x=ests)
    dmudtheta <- focus_deriv[1:pp]
    tau0sq <- t(dmudtheta) %*% solve(J00) %*% dmudtheta
    dmudgamma <- focus_deriv[pp + seq_len(qq)]
    omega <- J10 %*% solve(J00) %*% dmudtheta - dmudgamma
    psi.full <- t(omega) %*% deltahat

    if (sum(inds) > 0) { 
        Id <-  diag(rep(1,qq))
        Qinv <- solve(Q)
        pi.S <- matrix(Id[inds*(seq_len(qq)),],ncol=qq)

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
        var.adj.S <- tau0sq  +  t(omega) %*% Q0.S %*% omega 
    } else { 
        ## Special case for null model with all extra parameters excluded
        bias.S <- bias.adj.S <- psi.full
        var.adj.S <- tau0sq
        FIC.S <- bias.S^2
    }    
    mse.adj.S <- var.adj.S + bias.adj.S^2

    res <- c(FIC = FIC.S,
             rmse = sqrt(mse.adj.S / n),   # book p157
             bias = bias.S / sqrt(n),
             bias.adj = bias.adj.S / sqrt(n),
             se.adj = sqrt(var.adj.S / n)
             )
    res
}

FIC <- fic
