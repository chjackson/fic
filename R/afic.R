afic <- function(par, J, inds, inds0, gamma0=0, n, focus_deriv, wt)
{
    npar <- length(par)
    inds <- check_indsinds0(inds, inds0)
    i0 <- which(inds0==1)
    dnhat <- par[-i0] - gamma0
    J00 <- J[i0, i0, drop=FALSE]
    J10 <- J[-i0, i0, drop=FALSE]
    J01 <- J[i0,-i0, drop=FALSE]
    J11 <- J[-i0,-i0, drop=FALSE]
    Q <- solve(J)[-i0,-i0, drop=FALSE]   # using book notation.  called K in original code.
    dmudtheta <- focus_deriv[i0,,drop=FALSE]
    dmudgamma <- focus_deriv[-i0,,drop=FALSE]

    ## afic specific, could move all below
    ncov <- ncol(focus_deriv)
    BX <- array(dim=c(npar, npar, ncov))
    for (i in 1:ncov){
        dmu <- numeric(npar)
        dmu[i0] <- dmudtheta[,i]
        dmu[-i0] <- dmudgamma[,i]
        BX[,,i] <- outer(dmu, dmu)
    }
    B <- apply(BX, c(1,2), function(x)sum(x * wt))
    B00 <- B[i0, i0, drop=FALSE]
    B10 <- B[-i0, i0, drop=FALSE]
    B01 <- B[i0,-i0, drop=FALSE]
    B11 <- B[-i0,-i0, drop=FALSE]
    invJ00 <- solve(J00)
    A <- J10 %*% invJ00 %*% B00  %*% invJ00 %*% J01 -
      J10 %*% invJ00 %*% B01 -
      B10 %*% invJ00 %*% J01 +
      B11 
    tau0sqX <- diag(t(dmudtheta) %*% solve(J00) %*% dmudtheta)
    tau0sqA <- sum(tau0sqX*wt) # integrated over covariate space
    ## end afic specific 
    omega <- J10 %*% solve(J00) %*% dmudtheta - dmudgamma # q x m, where m is number of alternative focuses (typically covariate values)
    omegaA <- omega %*% wt

    indsS <- inds[inds0==0]
    qq <- sum(inds0==0)  # maximum number of "extra" parameters
    Id <-  diag(rep(1,qq))
    if (sum(indsS) > 0) {
        Qinv <- solve(Q)  # q x q 
        pi.S <- matrix(Id[indsS*(seq_len(qq)),],ncol=qq) # qs x q
        Q.S <- solve(pi.S %*% Qinv %*% t(pi.S)) # qs x qs
        Q0.S <- t(pi.S) %*% Q.S %*% pi.S        # q x q
        G.S <- Q0.S %*% Qinv # called M.S in original code.  q x q
        ### end 

        ### start afic specific 
        bias.S <- t(omegaA) %*% (Id - G.S) %*% dnhat
        sqbias <- sum(diag((Id - G.S) %*% outer(dnhat, dnhat) %*% t(Id - G.S) %*% A))
        IS <- sum(diag((Id - G.S) %*% (outer(dnhat, dnhat) - Q) %*% t(Id - G.S) %*% A))
        IIS <- sum(diag(Q0.S %*% A))  # diag(t(omega) %*% Q0.S %*% omega)
    } else {
        bias.S <- bias.adj.S <- t(omegaA) %*% dnhat
        sqbias <- bias.S^2
        IS <- sum(diag(Id %*% (outer(dnhat, dnhat) - Q) %*% t(Id) %*% A))
        IIS <- 0
    }

    var.S <- tau0sqA + IIS
    bias.adj.S <- sign(bias.S) * sqrt(max(IS, 0))
    afic_book <- bias.adj.S^2  + IIS # (6.27) in book

#    AFIC.S <- n*(sqbias + 2*IIS)   # quantity that reduces to standard FIC
#    mse.S <- sqbias + 2*IIS + tau0sqA - diag(t(omegaA) %*% Q %*% omegaA)
    mse.S <- IS + var.S
    AFIC.S <- n*(mse.S - tau0sqA + diag(t(omegaA) %*% Q %*% omegaA))
    mse.adj.S <- bias.adj.S^2 + var.S
    
    resA <- cbind(
        rmse     = sqrt_nowarning(mse.S),
        rmse.adj = sqrt_nowarning(mse.adj.S),   # book p157
        bias     = bias.adj.S,
        se       = sqrt(var.S),
        FIC      = AFIC.S
    )
    colnames(resA) <- c("FIC", "rmse", "rmse.adj", "bias", "se")

    resA
}


## We might be able to share more code between fic_core and AFIC,
## likewise for Cox models, but doesn't seem worth it unless we
## implement several more FIC-like methods.

