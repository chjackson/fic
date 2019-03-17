## same example as in vignette

context("Normal example")

wide.lm <- lm(mpg ~ am + wt + qsec + disp + hp, data=mtcars)
inds0 <- c(1, rep(0, length(coef(wide.lm))-1))
combs <- all_inds(wide.lm, inds0)

cmeans <- colMeans(model.frame(wide.lm)[,c("wt","qsec","disp","hp")])
X <- c(intercept=1, am=0, cmeans)
ficall <- fic(wide.lm, inds=combs, inds0=inds0, focus="mean_normal", X=X)
if (interactive()) ggplot_fic(ficall)
summary(ficall)
ficall[ficall$rmse <1,] # all include AM

wide.lm <- lm(mpg ~ cyl + disp + hp + drat + wt + qsec + vs + am + gear + carb, data=mtcars)

mod1.lm <- lm(mpg ~ cyl + disp + hp + drat + wt, data=mtcars)
inds0 <- c(1,0,0,0,0, 0,0,0,0,0,0)
inds1 <- c(1,1,1,1,1, 0,0,0,0,0,0)
X <- model.matrix(wide.lm)[1:5,]

ficall <- fic(wide.lm, inds=inds1, inds0=inds0, focus="mean_normal", X=X)
1

par <- coef(wide.lm)
J <- solve(vcov(wide.lm))
ana <- fic_multi(par=par, J=J,  inds=inds1, inds0=inds0, n=nrow(mtcars), focus="mean_normal", X=X, parsub=c(coef(mod1.lm), 0, 0, 0, 0, 0))
num <- fic_multi(par=par, J=J,  inds=inds1, inds0=inds0, n=nrow(mtcars), focus="mean_normal", X=X, parsub=c(coef(mod1.lm), 0, 0, 0, 0, 0))

test_that("Results of higher level and lower level functions match", {
    expect_equivalent(ana[,"FIC",1], ficall$FIC)
})
          
test_that("Analytic focus derivatives match numeric: normal", {
    expect_equivalent(ana, num)
})

combs <- all_inds(wide.lm, inds0)

X <- model.matrix(wide.lm)[1,]
ficres <- fic(wide.lm, combs, focus=mean_normal, X=X)
if (interactive()) ggplot_fic(ficres)


## Extra args to focus function 

focus_med <- function(par,X,sigma){
    qnorm(0.5, mean = as.numeric(X %*% par), sd=sigma)
}

focus_quantile <- function(par,X,sigma,focus_p=0.5){
    qnorm(focus_p, mean = as.numeric(X %*% par), sd=sigma)
}

test_that("focus functions with sigma and/or extra arguments",{
    ficres_med <- fic(wide.lm, combs[1:4,], focus=focus_med, X=X)
    expect_equal(ficres_med$FIC[1:4], ficres$FIC[1:4])
    ficres_q <- fic(wide.lm, combs[1:4,], focus=focus_quantile, X=X, focus_p=0.5)
    expect_equal(ficres_med$FIC[1:4], ficres_q$FIC[1:4])
})

test_that("focus function with multiple quantiles",{
    ficres_qmulti <- fic(wide.lm, combs[1:4,], focus=focus_quantile, X=X, focus_p=c(0.1, 0.5, 0.9))
    expect_equal(ficres$FIC[1:4], ficres_qmulti[ficres_qmulti$vals=="0.5","FIC"])
})

test_that("focus function argument length mismatch",{
    Xmat <- matrix(X, nrow=1)[c(1,1),]
    expect_error(fic(wide.lm, combs[1:4,], focus=focus_quantile, X=Xmat, focus_p=c(0.1, 0.5, 0.9)),
                 "Number of focuses")
})
