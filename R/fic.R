##' Focused Information Criterion: general function
##'
##' Focused information criterion for general models.
##'
##' @param ests Vector of maximum likelihood estimates from the wide model
##'
##' @param J Information matrix from the wide model, evaluated at the maximum likelihood estimates and divided by \code{n}.
##'
##' @param inds Vector of 0s and 1s of length \code{length(ests) - pp}, with 1s in the positions which define the submodel to be compared against the wide model.
##'
##' @param pp Number of parameters which we would always include in any submodel.  These parameters are assumed to be the first \code{pp} parameters included in \code{ests}.
##'
##' @param n Number of observations in the data used to fit the wide model.
##'
##' @param focus An R function taking one argument: a vector of parameters; and returning one scalar: the focus quantity of interest.  Not required if \code{focus.deriv} is specified
##'
##' @param focus.deriv Vector of partial derivatives of the focus function with respect to the parameters in the wide model.  If not supplied, this is computed from the function supplied in \code{focus}, using numerical differentiation.
##'
##' @return A vector containing the focused information criterion, bias and variance quantities (TO DEFINE) for the defined submodel.
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
                pp, #  number of covariates which we always include, ie length of theta (including intercept)
                n, #  observations 
                focus=NULL,
                focus.deriv=NULL
                ) #  function of wide params, returning focus 
{
    ## TODO: implement special case for sum(inds) = 0 
    deltahat <- sqrt(n)*ests[-(1:pp)]
    J00 <- J[1:pp, 1:pp]
    J10 <- J[-(1:pp), 1:pp]
    J01 <- J[1:pp,-(1:pp)]
    J11 <- J[-(1:pp),-(1:pp)]
    invJ <- solve(J)
    K <- invJ[-(1:pp),-(1:pp)]   # called Q in book.
    qq <- length(inds) # maximum number of "extra" covariates
    
    if(is.null(focus.deriv))
        focus.deriv <- numDeriv::grad(func=focus, x=ests)
    dmudtheta <- focus.deriv[1:pp]
    dmudgamma <- focus.deriv[pp + 1:qq]
    omega <- J10 %*% solve(J00) %*% dmudtheta - dmudgamma
    psi.full <- t(omega) %*% deltahat

    tau0sq = t(dmudtheta) %*% solve(J00) %*% dmudtheta

    Id <-  diag(rep(1,qq))
    Kinv <- solve(K)
    pi.S <- matrix(Id[inds*(1:qq),],ncol=qq)
    K.S <- solve(pi.S %*% Kinv %*% t(pi.S)) # called Q.S in book
    M.S <- t(pi.S) %*% K.S %*% pi.S %*% Kinv # called G.S in book
    psi.S <- t(omega)%*% M.S %*% deltahat
    omega.S <- pi.S %*% omega
    FIC.S <- (psi.full - psi.S)^2  +  2*t(omega.S) %*% K.S %*% omega.S

    bias.S <- t(omega) %*% deltahat - psi.S
    bias2.S <- sign(bias.S) * sqrt(max(0,t(omega) %*%
                                         (Id - M.S) %*%
                                         (deltahat %*% t(deltahat) - K) %*%
                                         (Id - M.S) %*%
                                         omega)) # sign*root of sqb3 on p152 
    ## V(S) on p152. or is it?  uses QO.S in book  = pis^t Q pis 
    var2.S <- tau0sq  +  t(omega) %*% M.S %*% omega    
    mse.S <- var2.S + bias2.S^2

    res <- c(FIC=FIC.S,
             bias=bias.S/n,
             bias2=bias2.S/n,
             var2=var2.S/n)
    res
}
