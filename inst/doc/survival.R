## ------------------------------------------------------------------------
if (!require("flexsurv")) stop("The `flexsurv` package should be installed to run code in this vignette") 
ex <-  flexsurvreg(Surv(recyrs, censrec) ~ 1, data=bc, dist="exponential")
we <- flexsurvreg(Surv(recyrs, censrec) ~ 1, data=bc, dist="weibull")
gg <- flexsurvreg(Surv(recyrs, censrec) ~ 1, data=bc, dist="gengamma")

## ------------------------------------------------------------------------
we2 <- flexsurvreg(Surv(recyrs, censrec) ~ 1, data=bc, dist="gengamma", inits=c(1,1,1), fixedpars=3)  # weibull with shape 1/sigma, scale exp(mu)
we3 <- flexsurvreg(Surv(recyrs, censrec) ~ 1, data=bc, dist="gengamma", inits=c(1,1,1), fixedpars=c(2,3)) # exponential model with rate 1/exp(mu)

## ------------------------------------------------------------------------
library(fic)
indmat <- rbind(ggamma = c(1,1,1),
                weib   = c(1,1,0),
                exp    = c(1,0,0))
gamma0 <- c(0,1)
focus <- function(par){
    rmst_gengamma(8, par[1], exp(par[2]), par[3])
}
fic(gg, inds=indmat, inds0=indmat[3,], gamma0=gamma0, focus=focus, sub=list(gg, we2, we3))

## ------------------------------------------------------------------------
focus(coef(gg))
focus(coef(we2))
focus(coef(we3))

## ------------------------------------------------------------------------
set.seed(1)
y <- rexp(100); cen <- rep(1,100)
gge <- flexsurvreg(Surv(y, cen) ~ 1, dist="gengamma")
focus <- function(par){
    pgengamma(1, par[1], exp(par[2]), par[3], lower.tail=FALSE)
}
gge
focus(coef(gge))
rbind(
    ggamma  = fic(gge, inds=indmat[1,], inds0=indmat[3,], gamma0=gamma0, focus=focus),
    weibull = fic(gge, inds=indmat[2,], inds0=indmat[3,], gamma0=gamma0, focus=focus),
    exp     = fic(gge, inds=indmat[3,], inds0=indmat[3,], gamma0=gamma0, focus=focus)
)


