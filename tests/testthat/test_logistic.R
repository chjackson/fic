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

test_that("Basic FIC and rmse results",{
    expect_equal(ficall$FIC, c(1.18660265602458, 1.30514383073585, 1.00470197274807), tol=1e-06)
    expect_equal(ficall$rmse, c(0.0722883108076118, 0.0803563210789494, 0.076428849974211), tol=1e-06)
    expect_equal(ficall$rmse.adj, c(0.0722883108076118, 0.0803563210789494, 0.076428849974211))
    expect_equal(ficall$bias, c(0.0459169526618119, 0.0730792045208761, 0.0610284223749518))
    expect_equal(ficall$se, c(0.0558321890818443, 0.0334150894647788, 0.0460097899452381))
    expect_equal(ficall$focus, c(0.397835568632103, 0.242744427996539, 0.320289998314321))
})

par <- coef(wide.glm)
n <- nrow(birthwt)
J <- solve(vcov(wide.glm))
ana <- fic_multi(par=par, J=J,  inds=inds1, inds0=inds0, n=n, focus="prob_logistic", X=X, parsub=c(coef(mod1.glm), 0, 0, 0, 0))
ana <- fic_multi(par=par, J=J,  inds=inds1, inds0=inds0, n=n, focus=prob_logistic, X=X, parsub=c(coef(mod1.glm), 0, 0, 0, 0))
num <- fic_multi(par=par, J=J,  inds=inds1, inds0=inds0, n=n, focus=focus_plogis, X=X, parsub=c(coef(mod1.glm), 0, 0, 0, 0))

test_that("Results of higher level and lower level functions match", {
    expect_equivalent(ana[,"FIC",1], ficall$FIC)
})
          
test_that("Analytic focus derivatives match numeric: normal", {
    expect_equivalent(ana, num)
})

test_that("Results of low and super-low level functions match",{
    focus_deriv <- prob_logistic_deriv(par=par, X=X)
    anacore <- fic_core(par=par, J=J, inds=inds1, inds0=inds0, gamma0=0, n=n, focus_deriv=focus_deriv)
    expect_equal(anacore[,"FIC"], ana[1:2,"FIC",1])
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

test_that("Fitting submodels",{
    ficall <- fic(wide.glm, inds=inds, inds0=inds0, focus=focus_plogis, X=X, sub=list(mod1.glm, mod2.glm))
    ficall2 <- fic(wide.glm, inds=inds, inds0=inds0, focus=focus_plogis, X=X)
    expect_equal(ficall$focus, ficall2$focus)
})

test_that("Supplying submodel parameters with multiple submodels",{
    expect_equivalent(focus_plogis(parsub[1,], X[2,]),
                 ficall$focus[ficall$mods=="inds1" & ficall$vals=="vals.nonsmoke"])
    expect_equivalent(focus_plogis(parsub[1,], X[2,]),
                      num["vals.nonsmoke", "focus", "inds1"])
})

#test_that("Narrow model parameters in the middle",{
inds0mid <- c(1,0,1,0,0,0,0,0)
fic(wide.glm, inds=inds1, inds0=inds0mid, focus=prob_logistic, X=X)
#})

test_that("Vector or matrix focuses allowed",{
    ficmat <- fic(wide.glm, inds=inds1, inds0=inds0mid, focus=prob_logistic, X=X)
    ficsm <- fic(wide.glm, inds=inds1, inds0=inds0mid, focus=prob_logistic, X=vals.smoke)
    ficnsm <- fic(wide.glm, inds=inds1, inds0=inds0mid, focus=prob_logistic, X=vals.nonsmoke)
    expect_equal(ficmat[ficmat$vals=="vals.smoke",c("rmse","focus")],
                 ficsm[,c("rmse","focus")])
})

test_that("Covariate weights", { 
    ficwt <- fic(wide.glm, inds=inds1, inds0=inds0, focus=prob_logistic, X=X, wt=c(0.1, 0.9))
    ficwt2 <- fic(wide.glm, inds=inds1, inds0=inds0, focus=prob_logistic, X=X, wt=c(0.9, 0.1))
    expect_lt(ficwt2$rmse[ficwt2$vals=="ave"], ficwt$rmse[ficwt$vals=="ave"])
})

test_that("Tidy and untidy output",{ 
    ficall <- fic(wide.glm, inds=inds1, inds0=inds0, focus=focus_plogis, X=X)
    ficall_untidy <- fic(wide.glm, inds=inds1, inds0=inds0, focus=focus_plogis, X=X, tidy=FALSE)
    expect_equivalent(ficall$rmse[ficall$vals=="vals.smoke"],
                      ficall_untidy["vals.smoke","rmse",])
})
