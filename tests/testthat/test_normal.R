## Normal example: focus is expected MPG for cars like the first two in the mtcars dataset

context("Normal example")

focus_mean <- function(par, X){X %*% par}
wide.lm <- lm(mpg ~ cyl + disp + hp + drat + wt + qsec + vs + am + gear + carb, data=mtcars)
mod1.lm <- lm(mpg ~ cyl + disp + hp + drat, data=mtcars)
inds0 <- c(1,0,0,0,0,0,0,0,0,0,0)
inds1 <- c(1,1,1,1,1,0,0,0,0,0,0)
X <- model.matrix(wide.lm)[1:5,]

ficall <- fic(wide.lm, inds=inds1, inds0=inds0, focus=focus_mean, X=X)

par <- coef(wide.lm)
J <- solve(vcov(wide.lm))/nrow(mtcars)
ana <- fic_multi(par=par, J=J,  inds=inds1, inds0=inds0, n=nrow(mtcars), focus="mean_normal", X=X, parsub=c(coef(mod1.lm), 0, 0, 0, 0, 0, 0))
num <- fic_multi(par=par, J=J,  inds=inds1, inds0=inds0, n=nrow(mtcars), focus=focus_mean, X=X, parsub=c(coef(mod1.lm), 0, 0, 0, 0, 0, 0))

test_that("Results of higher level and lower level functions match", {
    expect_equivalent(ana[,"FIC",1], ficall$FIC)
})
          
test_that("Analytic focus derivatives match numeric: normal", {
    expect_equivalent(ana, num)
})
