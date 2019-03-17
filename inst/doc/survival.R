## ------------------------------------------------------------------------
if (!require("flexsurv")) 
    stop("The `flexsurv` package should be installed
to run code in this vignette") 
ex <-  flexsurvreg(Surv(recyrs, censrec) ~ 1, data=bc, dist="exponential")
we <- flexsurvreg(Surv(recyrs, censrec) ~ 1, data=bc, dist="weibull")
gg <- flexsurvreg(Surv(recyrs, censrec) ~ 1, data=bc, dist="gengamma")

## ----fig.width=6,fig.height=5--------------------------------------------
plot(gg, ci=FALSE, conf.int=FALSE, ylab="Survival", xlab="Years")
lines(we, col="blue", ci=FALSE)
lines(ex, col="green", ci=FALSE)
legend("topright", lty=c(1,1,1), lwd=c(2,2,2), 
       col=c("green", "blue","red"), 
       c("Exponential","Weibull","Generalized gamma"))

## ------------------------------------------------------------------------
ex2 <- flexsurvreg(Surv(recyrs, censrec) ~ 1, data=bc, dist="gengamma", 
                   inits=c(1,1,1), fixedpars=c(2,3))

## ------------------------------------------------------------------------
we2 <- flexsurvreg(Surv(recyrs, censrec) ~ 1, data=bc, dist="gengamma", 
                   inits=c(1,1,1), fixedpars=3) 

## ----eval=FALSE----------------------------------------------------------
#  1 / exp(ex2$res["mu","est"])
#  1 / we2$res["sigma","est"]
#  exp(we2$est["mu","est"])

## ------------------------------------------------------------------------
focus <- function(par){
    rmst_gengamma(8, par[1], exp(par[2]), par[3])
}

## ------------------------------------------------------------------------
indmat <- rbind(exp    = c(1,0,0),
                weib   = c(1,1,0),
                ggamma = c(1,1,1))

## ------------------------------------------------------------------------
gamma0 <- c(0,1)

## ------------------------------------------------------------------------
library(fic)
fic(gg, inds=indmat, gamma0=gamma0, focus=focus, sub=list(ex2, we2, gg))

## ------------------------------------------------------------------------
set.seed(1)
y <- rexp(50); cen <- rep(1,50)
gge <- flexsurvreg(Surv(y, cen) ~ 1, dist="gengamma")
fic(gge, inds=indmat, gamma0=gamma0, focus=focus)

## ------------------------------------------------------------------------
if (!require("survival"))
    stop("The `survival` package should be installed
to run code in this vignette") 
ex <-  survreg(Surv(recyrs, censrec) ~ 1, data=bc, dist="exponential")
we <-  survreg(Surv(recyrs, censrec) ~ 1, data=bc, dist="weibull")
indmat <- rbind(exp    = c(1,0),
                weib   = c(1,1))

## ------------------------------------------------------------------------
focus <- function(par){
    rmst_weibull(8, 1/exp(par[2]), exp(par[1]))
}

## ------------------------------------------------------------------------
fic(we, inds=indmat, focus=focus, sub=list(ex, we))

