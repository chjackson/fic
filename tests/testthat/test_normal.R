## Normal example: focus is expected MPG for cars like the first two in the mtcars dataset

# MPG on TRANSMISSION TYPE, WEIGHT, and QUARTER MILE TIME 'may be best for prediction'. This model is difficult to interpret, and the absence of DISPLACEMENT or HORSEPOWER, which intuition suggests should be important in the prediction of MPG, is surprising. Hocking noted (pp. 25-26) that DISPLACEMENT exhibits instability and he found it disturbing that the ridge trace did not suggest the important role of his optimal subset

context("Normal example")

wide.lm <- lm(mpg ~ am + wt + qsec + disp + hp, data=mtcars)
inds0 <- c(1, rep(0, length(coef(wide.lm))-1))
combs <- all_inds(wide.lm, inds0)

cmeans <- colMeans(model.frame(wide.lm)[,c("wt","qsec","disp","hp")])
X <- c(intercept=1, am=0, cmeans)
ficall <- fic(wide.lm, inds=combs, inds0=inds0, focus="mean_normal", X=X)
ggplot_fic(ficall)
summary(ficall)
ficall[ficall$rmse <1,] # all include AM

## todo plot against AIC / BIC 

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
(ficres <- fic(wide.lm, combs, focus=mean_normal, X=X))
ggplot_fic(ficres)

require(graphics)
pairs(mtcars, main = "mtcars data")
coplot(mpg ~ disp | as.factor(cyl), data = mtcars,
       panel = panel.smooth, rows = 1)

## TODO gapminder

## extra args to focus function 

focus_med <- function(par,X,sigma){
    qnorm(0.5, mean = as.numeric(par %*% X), sd=sigma)
}

focus_quantile <- function(par,X,sigma,focus_p=0.5){
    qnorm(focus_p, mean = as.numeric(par %*% X), sd=sigma)
}

## expect that fic same for mean and any quantile 

