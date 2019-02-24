context("Error handling")

wide.glm <- glm(low ~ lwtkg + age + smoke + ht + ui + smokeage + smokeui, data=birthwt, family=binomial)
mod1.glm <- glm(low ~ lwtkg + age + smoke, data=birthwt, family=binomial)
inds0 <- c(1,1,0,0,0,0,0,0)
inds1 <- c(1,1,1,1,0,0,0,0)
focus_plogis <- function(par, X)plogis(X %*% par)
vals.smoke <-    c(1, 58.24, 22.95, 1, 0, 0, 22.95, 0)
vals.nonsmoke <- c(1, 59.50, 23.43, 0, 0, 0, 0, 0)
X <- rbind(vals.smoke, vals.nonsmoke)


expect_error(fic(wide=wide.glm, inds=c(inds1,0), inds0=inds0, focus=focus_plogis, X=X),
             "`inds` of length")
expect_error(fic(wide=wide.glm, inds=inds1, inds0=inds0[-1], focus=focus_plogis, X=X),
             "Length of `inds0` must match number of parameters")
expect_error(fic(wide.glm, inds=inds1, inds0=inds0, focus="foo", X=X), "focus function `foo` not found")

Xwrong <- cbind(X, c(1,1))
expect_error(fic(wide.glm, inds=inds1, inds0=inds0, focus=focus_plogis, X=Xwrong),
             "Number of columns of X is")

expect_error(fic(wide.glm, inds=inds1, inds0=inds0, focus=focus_plogis, X=X, gamma0=rep(0, 7)),
             "Length of gamma0")

expect_error(
    fic(2, inds=inds1, inds0=inds0, focus=focus_plogis, X=X, gamma0=rep(0, 7)),
    "argument is specified correctly") ## TODO nicer way to detect that something is a fitted model object.  par, n and J should be extracted with no error and be of right form


expect_error(fic(wide.glm, inds=inds1, inds0=inds0, focus=focus_plogis, X=X, sub=mod1.glm),
             "`sub` should be a list of fitted model objects")

expect_error(fic(wide.glm, inds=inds1, inds0=inds0, focus=focus_plogis, X=X, sub=list(mod1.glm, mod2.glm)), "not found")
mod2.glm <- glm(low ~ lwtkg + age, data=birthwt, family=binomial)
expect_error(fic(wide.glm, inds=inds1, inds0=inds0, focus=focus_plogis, X=X, sub=list(mod1.glm, mod2.glm)), "`sub` of length")


## TODO errors for parsub in fic_multi

## Lengths of parameters in submodels don't match inds

## errors in lower-level functions 

## TODO mismatch between X and length of par
## Currently check_X checks for X too big
## Hard to check for X too small, as covariate effects might be a subset of full set of pars, and hard to identify which pars are the covariate effects 
## This check might need to be done within the focus function. 

expect_error(get_fns(list(foo=1, bar=2)), "components named")
expect_error(get_fns(list(coef=1, bar=2)), "components named")
expect_error(get_fns(list(coef=1)), "should be a function")



## TODO fns argument not supplied for new model class 
