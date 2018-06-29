afic <- function(par, J, inds, inds0, gamma0=0, n, focus_deriv, Xwt)
{
    npar <- length(par)
    inds <- check_indsinds0(inds, inds0)
    ncov <- ncol(focus_deriv)
    i0 <- which(inds0==1)
    deltahat <- sqrt(n)*(par[-i0] - gamma0)
    J00 <- J[i0, i0, drop=FALSE]
    J10 <- J[-i0, i0, drop=FALSE]
    J01 <- J[i0,-i0, drop=FALSE]
    J11 <- J[-i0,-i0, drop=FALSE]
    Q <- solve(J)[-i0,-i0, drop=FALSE]   # using book notation.  called K in original code.
    dmudtheta <- focus_deriv[i0,,drop=FALSE]
    dmudgamma <- focus_deriv[-i0,,drop=FALSE]
    BX <- array(dim=c(npar, npar, ncov))
    for (i in 1:ncov){
        dmu <- c(dmudtheta[,i], dmudgamma[,i])
        BX[,,i] <- outer(dmu, dmu)
    }

    B <- apply(BX, c(1,2), function(x)sum(x * Xwt))
    B00 <- B[i0, i0, drop=FALSE]
    B10 <- B[-i0, i0, drop=FALSE]
    B01 <- B[i0,-i0, drop=FALSE]
    B11 <- B[-i0,-i0, drop=FALSE]    
    invJ00 <- solve(J00)
    A <- J10 %*% invJ00 %*% B00  %*% invJ00 %*% J01 -
      J10 %*% invJ00 %*% B01 -
      B10 %*% invJ00 %*% J01 +
      B11 

    indsS <- inds[inds0==0]
    qq <- sum(inds0==0)  # maximum number of "extra" parameters
    Id <-  diag(rep(1,qq))
    tau0sqX <- diag(t(dmudtheta) %*% solve(J00) %*% dmudtheta)
    tau0sq <- sum(tau0sqX*Xwt) # integrated over covariate space
    omega <- J10 %*% solve(J00) %*% dmudtheta - dmudgamma # q x m, where m is number of alternative focuses (typically covariate values)
    omega <- omega %*% Xwt
    if (sum(indsS) > 0) {
        Qinv <- solve(Q)  # q x q 
        pi.S <- matrix(Id[indsS*(seq_len(qq)),],ncol=qq) # qs x q
        Q.S <- solve(pi.S %*% Qinv %*% t(pi.S)) # qs x qs
        Q0.S <- t(pi.S) %*% Q.S %*% pi.S        # q x q
        G.S <- Q0.S %*% Qinv # called M.S in original code.  q x q
        bias.S <- t(omega) %*% (Id - G.S) %*% deltahat
        sqbias <- sum(diag((Id - G.S) %*% outer(deltahat, deltahat) %*% t(Id - G.S) %*% A))
        IS <- sum(diag((Id - G.S) %*% (outer(deltahat, deltahat) - Q) %*% t(Id - G.S) %*% A))
        IIS <- sum(diag(Q0.S %*% A))  # diag(t(omega) %*% Q0.S %*% omega)
    } else {
        bias.S <- bias.adj.S <- t(omega) %*% deltahat
        sqbias <- bias.S^2
        var.S <- tau0sq
        IS <- sum(diag(Id %*% (outer(deltahat, deltahat) - Q) %*% t(Id) %*% A))
        IIS <- 0
    }
    var.S <- tau0sq + IIS
    bias.adj.S <- sqrt(max(IS, 0))
    afic_book <- bias.adj.S^2  + IIS # (6.27) in book
    AFIC.S <- sqbias + 2*IIS   # quantity that reduces to standard FIC
    mse.S <- AFIC.S + tau0sq - diag(t(omega) %*% Q %*% omega)
    mse.adj.S <- bias.adj.S^2 + var.S
    
    res <- cbind(
        FIC      = AFIC.S,
        rmse     = sqrt(mse.S / n),
        rmse.adj = sqrt(mse.adj.S / n),   # book p157
        bias     = bias.S / sqrt(n),
        bias.adj = bias.adj.S / sqrt(n),
        se       = sqrt(var.S / n)
    )
    colnames(res) <- c("FIC", "rmse", "rmse.adj", "bias", "bias.adj", "se")
    res
}


## TODO can we share more code between fic_core and AFIC 
