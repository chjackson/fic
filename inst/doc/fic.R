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
prob_logistic <- function(par, X)plogis(X %*% par)
vals.smoke <-    c(1, 58.24, 22.95, 1, 0, 0, 22.95, 0)
vals.nonsmoke <- c(1, 59.50, 23.43, 0, 0, 0, 0, 0)
X <- rbind("Smokers"=vals.smoke, "Non-smokers"=vals.nonsmoke)

## ---------------------------------------------------------
prob_logistic(coef(wide.glm), X=X)

## ---------------------------------------------------------
mod1.glm <- glm(low ~ lwtkg + age + smoke, data=birthwt, family=binomial)
mod2.glm <- glm(low ~ lwtkg + age + smoke + ht, data=birthwt, family=binomial)

## ---------------------------------------------------------
inds <- rbind(mod1 = c(1,1,1,1,0,0,0,0),
              mod2 = c(1,1,1,1,1,0,0,0))
inds0 <- c(1,1,0,0,0,0,0,0)

## ---------------------------------------------------------
fic1 <- fic(wide=wide.glm, inds=inds, inds0=inds0, focus=prob_logistic, X=X)
fic1

## ---------------------------------------------------------
combs <- all_inds(wide.glm, inds0)

## ---------------------------------------------------------
combs <- with(combs,
              combs[!((smoke==0 & smokeage==1) |
                      (smoke==0 & smokeui==1) |
                      (age==0 & smokeage==1) |
                      (ui==0 & smokeui==1)),])
ficres <- fic(wide=wide.glm, inds=combs, inds0=inds0,
              focus=prob_logistic, X=X)

## ---------------------------------------------------------
ggplot_fic(ficres)

## ---------------------------------------------------------
summary(ficres)

## ----eval=FALSE-------------------------------------------
#  fns <- list(coef = function(x)coef(x),
#              nobs = function(x)nobs(x),
#              vcov = function(x)vcov(x))
#  fic1 <- fic(wide=wide.glm, inds=inds, inds0=inds0, focus=prob_logistic,
#              fns=fns, X=X)

## ---------------------------------------------------------
library(survival)
wide <- coxph(Surv(years, death==1) ~ sex + thick_centred + infilt + epith + 
                ulcer + depth + age, data=melanoma)

## ---------------------------------------------------------
inds0 <- expand_inds(c(1,0,0,0,0,0,0), wide)
inds0

## ----eval=FALSE-------------------------------------------
#  list(focus=fic:::cox_hr, deriv=fic:::cox_hr_deriv,
#       dH=fic:::cox_hr_dH)
#  list(focus=fic:::cox_cumhaz, deriv=fic:::cox_cumhaz_deriv,
#       dH=fic:::cox_cumhaz_dH)
#  list(focus=fic:::cox_survival, deriv=fic:::cox_survival_deriv,
#       dH=fic:::cox_survival_dH)

## ---------------------------------------------------------
combs <- all_inds(wide,inds0)

## ---------------------------------------------------------
newdata <- with(melanoma,
                data.frame(sex = c("female","male"),
                           thick_centred = tapply(thick_centred, sex, mean),
                           infilt=4, epith=1, ulcer=1, depth=2,
                           age = tapply(age, sex, mean)))
X <- newdata_to_X(newdata, wide, intercept=FALSE)
ficall <- fic(wide, inds=combs, inds0=inds0, focus="survival", X=X, t=5)
ggplot_fic(ficall, ci=FALSE, xlim=c(0,1))

## ---------------------------------------------------------
summary(ficall)

## ----echo=FALSE,eval=FALSE--------------------------------
#  par(mfrow=c(1,2))
#  plot(ficall$FIC[ficall$vals=="female"], ficall$focus[ficall$vals=="female"],
#       xlim=c(0,30), ylim=c(0, 1),
#       ylab = "5yr survival estimates, women", xlab="FIC")
#  plot(ficall$FIC[ficall$vals=="male"], ficall$focus[ficall$vals=="male"],
#       xlim=c(0,30), ylim=c(0.2, 0.9),
#       ylab = "5yr survival estimates, men", xlab="FIC")

