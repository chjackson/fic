## ------------------------------------------------------------------------

if (!require("sn")) stop("The `sn` package should be installed to run code in this vignette") 
data(ais)
plot(density(ais$Hc))
plot(ais$BMI, ais$Hc, pch=19)

ldsnorm <- function(x, mean, sd, skew){
    log(skew) + (skew-1)*pnorm(x, mean, sd, log=TRUE) + dnorm(x, mean, sd, log=TRUE)
}

mloglik <- function(b0, b1, sd, skew){
    ret <- -sum(ldsnorm(ais$Hc, b0 + b1*ais$BMI, sd, skew))
    ret
}

snlm <- function(reg=FALSE, skew=FALSE){
    inilm <- lm(Hc ~ BMI, data=ais)
    ini <- c(coef(inilm), log(summary(inilm)$sigma), 0)
    if (!reg && !skew){
        ini <- ini[c(1,3)]
        fn <- function(par) mloglik(par[1], 0, exp(par[2]), 1)
    } else if (reg && !skew){
        ini <- ini[c(1,2,3)]
        fn <- function(par) mloglik(par[1], par[2], exp(par[3]), 1)
    } else if (!reg && skew){
        ini <- ini[c(1,3,4)]
        fn <- function(par) mloglik(par[1], 0, exp(par[2]), exp(par[3]))
    } else if (reg && skew){
        fn <- function(par) mloglik(par[1], par[2], exp(par[3]), exp(par[4]))
    }
    opt <- nlm(fn, ini, hessian=TRUE)
    vcov <- solve(opt$hessian)
    list(loglik=-opt$minimum, est=opt$estimate, vcov=vcov, nobs=nrow(ais))
}


## ------------------------------------------------------------------------


(mod1 <- snlm(reg=FALSE, skew=FALSE))
(mod2 <- snlm(reg=TRUE, skew=FALSE))
(mod3 <- snlm(reg=FALSE, skew=TRUE))
(mod4 <- snlm(reg=TRUE, skew=TRUE))



## ------------------------------------------------------------------------

focus <- function(par, X){
    X %*% par[1:2]
}

fns <- list(coef=function(x)x$est,
            vcov=function(x)x$vcov,
            nobs=function(x)x$nobs)

med.bmi <- rbind(male=c(1, 23.56), female=c(1, 21.82))

inds <- rbind(c(1,0,1), c(1,1,1))
fic(mod2, inds=inds, inds0=c(1,0,1), fns=fns, focus=focus, X=med.bmi,
    sub=list(mod1, mod2))




## ------------------------------------------------------------------------

inds <- rbind(c(1,0,1,0), c(1,1,1,0), c(1,0,1,1), c(1,1,1,1))



