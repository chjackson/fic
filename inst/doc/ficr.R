## ------------------------------------------------------------------------
library(fic)
wide.glm <- glm(low ~ lwtkg + age + smoke + ht + ui + smokeage + smokeui, data=birthwt, family=binomial)

## ------------------------------------------------------------------------
vals.smoke <-    c(1, 58.24, 22.95, 1, 0, 0, 22.95, 0)
focus.smoke <- function(par)plogis(par %*% vals.smoke)

## ------------------------------------------------------------------------
vals.nonsmoke <- c(1, 59.50, 23.43, 0, 0, 0, 0, 0)
focus.nonsmoke <- function(par)plogis(par %*% vals.nonsmoke)

## ------------------------------------------------------------------------
focus.smoke(coef(wide.glm))
focus.nonsmoke(coef(wide.glm))

## ------------------------------------------------------------------------
mod1.glm <- glm(low ~ lwtkg + age + smoke, data=birthwt, family=binomial)
mod2.glm <- glm(low ~ lwtkg + age + smoke + ht, data=birthwt, family=binomial)

## ------------------------------------------------------------------------
inds0 <- c(1,1,0,0,0,0,0,0)
inds1 <- c(1,1,1,1,0,0,0,0)
inds2 <- c(1,1,1,1,1,0,0,0)
fic1 <- fic(wide=wide.glm, sub=mod1.glm, inds=inds1, inds0=inds0, focus=focus.smoke)
fic2 <- fic(wide=wide.glm, sub=mod2.glm, inds=inds2, inds0=inds0, focus=focus.smoke)
res <- rbind("No HT"=fic1, "HT"=fic2)
res

## ------------------------------------------------------------------------
par <- coef(wide.glm)
J <- solve(vcov(wide.glm))/nrow(birthwt)
fic(par=par, J=J,  inds=inds1, inds0=inds0, n=nrow(birthwt), focus=focus.smoke, parsub=coef(mod1.glm))

