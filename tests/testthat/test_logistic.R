context("Logistic regression")

wide.glm <- glm(low ~ lwtkg + age + smoke + ht + ui + smokeage + smokeui, data=birthwt, family=binomial)
mod1.glm <- glm(low ~ lwtkg + age + smoke, data=birthwt, family=binomial)
inds0 <- c(1,1,0,0,0,0,0,0)
inds1 <- c(1,1,1,1,0,0,0,0)
focus_plogis <- function(par, X)plogis(X %*% par)
vals.smoke <-    c(1, 58.24, 22.95, 1, 0, 0, 22.95, 0)
vals.nonsmoke <- c(1, 59.50, 23.43, 0, 0, 0, 0, 0)
X <- rbind(vals.smoke, vals.nonsmoke)

ficall <- fic(wide.glm, inds=inds1, inds0=inds0, focus=focus_plogis, X=X)

par <- coef(wide.glm)
n <- nrow(birthwt)
J <- solve(vcov(wide.glm))
ana <- fic_multi(par=par, J=J,  inds=inds1, inds0=inds0, n=n, focus="prob_logistic", X=X, parsub=c(coef(mod1.glm), 0, 0, 0, 0))
num <- fic_multi(par=par, J=J,  inds=inds1, inds0=inds0, n=n, focus=focus_plogis, X=X, parsub=c(coef(mod1.glm), 0, 0, 0, 0))

test_that("Results of higher level and lower level functions match", {
    expect_equivalent(ana[,"FIC",1], ficall$FIC)
})
          
test_that("Analytic focus derivatives match numeric: normal", {
    expect_equivalent(ana, num)
})

test_that("Supplying submodel parameters with single submodel",{
    ficall <- fic(wide.glm, inds=inds1, inds0=inds0, focus=focus_plogis, X=X, sub=list(mod1.glm))
    expect_equivalent(ficall$focus, num[,"focus",])
})

mod2.glm <- glm(low ~ lwtkg + age + smoke + ht, data=birthwt, family=binomial)
inds2 <- c(1,1,1,1,1,0,0,0)
inds <- rbind(inds1, inds2)

num <- fic_multi(par=par, J=J,  inds=inds, inds0=inds0, n=n, focus=focus_plogis, X=X)
ficall <- fic(wide.glm, inds=inds, inds0=inds0, focus=focus_plogis, X=X)

test_that("Multiple submodels",{
    expect_equivalent(as.vector(t(num[,"FIC",])), ficall$FIC)
})

parsub <- rbind(c(coef(mod1.glm), 0, 0, 0, 0),
                c(coef(mod2.glm),    0, 0, 0))
num <- fic_multi(par=par, J=J,  inds=inds, inds0=inds0, n=n, focus=focus_plogis, X=X, parsub=parsub)
ficall <- fic(wide.glm, inds=inds, inds0=inds0, focus=focus_plogis, X=X, sub=list(mod1.glm, mod2.glm))

test_that("Supplying submodel parameters with multiple submodels",{
    expect_equivalent(focus_plogis(parsub[1,], X[2,]),
                 ficall$focus[ficall$mods=="inds1" & ficall$vals=="vals.nonsmoke"])
    expect_equivalent(focus_plogis(parsub[1,], X[2,]),
                      num["vals.nonsmoke", "focus", "inds1"])
})



#### TESTS TODO

## Narrow model parameters in the middle


## sum(indsS) == 0

## calling core with multiple vector focuses gives same as matrix 
