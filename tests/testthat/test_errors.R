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
expect_error(fic(wide.glm, inds=inds1, inds0=inds0, focus="foo", X=X), "not found")

Xwrong <- cbind(X, c(1,1))
expect_error(fic(wide.glm, inds=inds1, inds0=inds0, focus=focus_plogis, X=Xwrong),
             "focus function returned an error")

expect_error(fic(wide.glm, inds=inds1, inds0=inds0, focus=focus_plogis, X=X, gamma0=rep(0, 7)),
             "Length of gamma0")

expect_error(
    fic(2, inds=inds1, inds0=inds0, focus=focus_plogis, X=X, gamma0=rep(0, 7)),
    "argument is specified correctly")


expect_error(fic(wide.glm, inds=inds1, inds0=inds0, focus=focus_plogis, X=X, sub=mod1.glm),
             "`sub` should be a list of fitted model objects")

expect_error(
    fic(wide.glm, inds=inds1, inds0=inds0, focus=focus_plogis, X=X, sub=list(mod1.glm, mod2.glm))
  , "not found")
mod2.glm <- glm(low ~ lwtkg + age, data=birthwt, family=binomial)
expect_error(fic(wide.glm, inds=inds1, inds0=inds0, focus=focus_plogis, X=X, sub=list(mod1.glm, mod2.glm)), "`sub` of length")

## Mismatch between X and length of par
Xbig <- X[,c(1:8,8)]
expect_error(
    fic(wide=wide.glm, inds=inds1, inds0=inds0, focus=focus_plogis, X=Xbig),
    "focus function returned an error")

expect_error(get_fns(list(foo=1, bar=2)), "components named")
expect_error(get_fns(list(coef=1, bar=2)), "components named")
expect_error(get_fns(list(coef=1)), "should be a function")

test_that("wt of wrong length",{
    expect_error(fic(wide.glm, inds=inds1, inds0=inds0, focus=prob_logistic, X=X, wt=c(0.1, 0.9, 0.1)), "of length 3")
    expect_error(fic(wide.glm, inds=inds1, inds0=inds0, focus=prob_logistic, X=X, wt=c(0.1)), "of length 1")
})

test_that("wrong argument names in focus",{
    focus <- function(ests){
        plogis(q = ests %*% vals.first)
    }
    expect_error(fic(wide.glm, inds=inds1, inds0=inds0, focus=focus, X=X),
                 "First argument of focus function")
})

test_that("args wrong way round in focus calculation",{
    focus <- function(par,X){
        plogis(q = par %*% X)
    }
    expect_error(
        fic(wide.glm, inds=inds1, inds0=inds0, focus=focus, X=X),
        "focus function returned an error")
})

## Error could be more specific, could check in check_focus that all args with no defaults have been supplied 
### Could check that any mandatory args with no defaults have been supplied. 
## give error if anames contains X, and has no default, but it's not in eargs 
## maybe not, because can't check for use of missing() 
# fic(wide.glm, inds=inds1, inds0=inds0, focus=focus)

pars <- coef(wide.glm)
n <- nrow(birthwt)
J <- solve(vcov(wide.glm))

test_that("parsub wrong length in lower level functions",{
expect_error( 
    fic_multi(par=pars, J=J,  inds=inds1, inds0=inds0, n=n, focus="prob_logistic", X=X,
              parsub=c(coef(mod1.glm), 0, 0, 0, 0, 0)),
    "parsub of length 9, should be 8")
expect_error( 
    fic_multi(par=pars, J=J,  inds=inds1, inds0=inds0, n=n, focus="prob_logistic", X=X,
              parsub=c(coef(mod1.glm), 0, 0, 0)),
    "parsub of length 7, should be 8")
})

test_that("wrong inds dimension in expand_inds", { 
    bwt.glm <- glm(low ~ lwtkg + age + smoke + ftv, data=birthwt, family="binomial")
    inds <- rbind(c(1,1,0,0),
                  c(1,1,1,1))
    expect_error(expand_inds(inds, bwt.glm),
                 "`inds` has 4 columns, but 5 terms in model")
    inds <- c(1,1,0,0)
    expect_error(expand_inds(inds, bwt.glm),
                 "`inds` of length 4, but 5 terms in model")
})

test_that("wrong inds dimension in all_inds", {
    bwt.glm <- glm(low ~ lwtkg + age + smoke, data=birthwt, family="binomial")
    expect_error(all_inds(bwt.glm, inds0=c(1,0,0)),
                 "inds0 of length 3, but 4 parameters in wide model")
})
