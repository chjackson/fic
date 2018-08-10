## ----echo=FALSE-------------------------------------------
options(width=60,digits=3)
options(prompt="R> ")
library(knitr)
opts_chunk$set(fig.path="fic-")

## ---------------------------------------------------------
library(fic)
wide.glm <- glm(low ~ lwtkg + age + smoke + ht + ui + smokeage + smokeui, 
                data=birthwt, family=binomial)

## ---------------------------------------------------------
focus <- function(par, X)plogis(X %*% par)
vals.smoke <-    c(1, 58.24, 22.95, 1, 0, 0, 22.95, 0)
vals.nonsmoke <- c(1, 59.50, 23.43, 0, 0, 0, 0, 0)
X <- rbind("Smokers"=vals.smoke, "Non-smokers"=vals.nonsmoke)

## ---------------------------------------------------------
focus(coef(wide.glm), X=X)

## ---------------------------------------------------------
mod1.glm <- glm(low ~ lwtkg + age + smoke, data=birthwt, family=binomial)
mod2.glm <- glm(low ~ lwtkg + age + smoke + ht, data=birthwt, family=binomial)

## ---------------------------------------------------------
inds <- rbind(mod1 = c(1,1,1,1,0,0,0,0),
              mod2 = c(1,1,1,1,1,0,0,0))
inds0 <- c(1,1,0,0,0,0,0,0)

## ---------------------------------------------------------
sub <- list(mod1.glm, mod2.glm)
fic1 <- fic(wide=wide.glm, inds=inds, inds0=inds0, focus=focus, X=X, sub=sub)
fic1

## ---------------------------------------------------------
combs <- all_inds(wide.glm, inds0)

## ---------------------------------------------------------
combs <- with(combs,
              combs[!((smoke==0 & smokeage==1) |
                      (smoke==0 & smokeui==1) |
                      (age==0 & smokeage==1) |
                      (ui==0 & smokeui==1)),])
ficres <- fic(wide=wide.glm, inds=combs, inds0=inds0, focus=focus, X=X)

## ----warning=FALSE----------------------------------------
nmod <- nrow(combs)
sub <- vector(nmod, mode="list")
XZ <- model.matrix(wide.glm)
for (i in 1:nmod){
  XZi <- XZ[,which(combs[i,]==1)]
  sub[[i]] <- glm(low ~ XZi - 1, data=birthwt, family=binomial)
}
ficres <- fic(wide=wide.glm, inds=combs, inds0=inds0, focus=focus, X=X, 
              sub=sub)

## ---------------------------------------------------------
ggplot_fic(ficres)

## ----eval=FALSE-------------------------------------------
#  fns <- list(coef = function(x)coef(x),
#              nobs = function(x)nobs(x),
#              vcov = function(x)vcov(x))
#  fic1 <- fic(wide=wide.glm, inds=inds, inds0=inds0, focus=focus, fns=fns,
#              X=X, sub=sub)

## ---------------------------------------------------------
wide <- coxph(Surv(years, death==1) ~ sex + thick_centred + infilt + epith + 
                ulcer + depth + age, data=melanoma)

## ---------------------------------------------------------
inds0 <- expand_inds(c(1,0,0,0,0,0,0), wide)
inds0

## ---------------------------------------------------------
combs <- all_inds(wide,inds0,intercept=FALSE)
nmod <- nrow(combs)
sub <- vector(nmod, mode="list")
for (i in 1:nmod){
  XZi <- model.matrix(wide)[,which(combs[i,]==1)]
  sub[[i]] <- coxph(Surv(years, death==1) ~ XZi - 1, data=melanoma)
}

## ---------------------------------------------------------
newdata <- with(melanoma,
                data.frame(sex = c("female","male"),
                           thick_centred = tapply(thick_centred, sex, mean),
                           infilt=4, epith=1, ulcer=1, depth=2,
                           age = tapply(age, sex, mean)))
X <- newdata_to_X(newdata, wide)
ficall <- fic(wide, inds=combs, inds0=inds0, focus="survival", X=X, t=5, sub=sub)
plot(ficall, xlim=c(0,1), ci=FALSE)
ggplot_fic(ficall, ci=FALSE)

