## ------------------------------------------------------------------------
ldsnorm <- function(x, mean, sd, lambda){
    log(lambda) + (lambda-1)*pnorm(x, mean, sd, log.p=TRUE) +
        dnorm(x, mean, sd, log=TRUE)
}

## ----fig.height=5--------------------------------------------------------
if (!require("sn"))
    stop("The `sn` package should be installed to run code in this vignette") 
data(ais)
par(mfrow=c(1,2))
plot(density(ais$Hc), xlab="Haematocrit level", main="")
plot(ais$BMI, ais$Hc, pch=19,
     xlab="Body mass index", ylab="Haematocrit level")

## ------------------------------------------------------------------------
mloglik <- function(b0, b1, sd, lambda){
    -sum(ldsnorm(ais$Hc, b0 + b1*ais$BMI, sd, lambda))
}

## ------------------------------------------------------------------------
fn1 <- function(par) mloglik(par[1], 0, exp(par[2]), 1)
fn2 <- function(par) mloglik(par[1], par[2], exp(par[3]), 1)
fn3 <- function(par) mloglik(par[1], 0, exp(par[2]), exp(par[3]))
fn4 <- function(par) mloglik(par[1], par[2], exp(par[3]), exp(par[4]))

## ------------------------------------------------------------------------
lm2 <- lm(Hc ~ BMI, data=ais)
cf <- unname(coef(lm2))
ini <- c(beta0=cf[1], beta1=cf[2], 
         logsigma=log(summary(lm2)$sigma), loglambda=0)

## ----warning=FALSE-------------------------------------------------------
opt1 <- nlm(fn1, ini[c("beta0","logsigma")], hessian=TRUE)
opt2 <- nlm(fn2, ini[c("beta0","beta1","logsigma")], hessian=TRUE)
opt3 <- nlm(fn3, ini[c("beta0","logsigma","loglambda")], hessian=TRUE)
opt4 <- nlm(fn4, ini, hessian=TRUE)

## ------------------------------------------------------------------------
mod1 <- list(est=opt1$estimate, vcov=solve(opt1$hessian) )
mod2 <- list(est=opt2$estimate, vcov=solve(opt2$hessian) )
mod3 <- list(est=opt3$estimate, vcov=solve(opt3$hessian) )
mod4 <- list(est=opt4$estimate, vcov=solve(opt4$hessian) )

## ------------------------------------------------------------------------
mean_snorm <- function(mu, sigma, lambda){
    f <- function(u){u*exp(ldsnorm(u, 0, 1, lambda))}
    mu + sigma * integrate(f, -Inf, Inf)$value
}
median_snorm <- function(mu, sigma, lambda){
    mu + sigma * qnorm(0.5^(1/lambda))
}

## ------------------------------------------------------------------------
focus1 <- function(par, X){
    mean_snorm(mu = par[1] + X %*% par[2],
               sigma = exp(par[3]), lambda = exp(par[4]))
}
focus2 <- function(par, X){
    median_snorm(mu = par[1] + X %*% par[2],
                 sigma = exp(par[3]), lambda = exp(par[4]))
}

## ------------------------------------------------------------------------
inds <- rbind("intcpt"     =c(1,0,1,0),
              "cov"        =c(1,1,1,0),
              "intcpt_skew"=c(1,0,1,1),
              "cov_skew"   =c(1,1,1,1))

fns <- list(coef=function(x)x$est,
            vcov=function(x)x$vcov,
            nobs=function(x)nrow(ais))

med.bmi <- rbind(male=23.56,  female=21.82)

library(fic)
fmean <- fic(mod4, inds=inds, fns=fns, focus=focus1, X=med.bmi, FIC=TRUE,
             sub=list(mod1, mod2, mod3, mod4))
fmean 
fmed <- fic(mod4, inds=inds, fns=fns, focus=focus2, X=med.bmi, FIC=TRUE,
            sub=list(mod1, mod2, mod3, mod4))
fmed

## ----fig.width=7,fig.height=4--------------------------------------------
ggplot_fic(fmean)
ggplot_fic(fmed)

