#library(devtools)
#load_all("../../msm/msm/")
#library(msm)
#load_all("fic")

if (requireNamespace("msm", quietly = TRUE)){

## 4-state PsA model using "psor" data in msm
## Wide model includes two covariates on all rates 
psor.q <- rbind(c(0,0.1,0,0), c(0,0,0.1,0), c(0,0,0,0.1), c(0,0,0,0))
psor.wide.msm <- msm::msm(state ~ months, subject=ptnum, data=psor, qmatrix = psor.q,  covariates = ~ollwsdrt+hieffusn, control=list(fnscale=1))
psor.wide.msm

## focus: total length of stay in state 4 up to 10 years
focus <- function(pars){
    x.new <- msm::updatepars.msm(psor.wide.msm, pars)
    msm::totlos.msm(x.new, covariates=0, t=10)["State 4"]
}

msm::totlos.msm(psor.wide.msm, covariates=0, t=10)
focus(psor.wide.msm$estimates) # should match 

## FIC analysis: omit covariate effects in turn from full model 
## Why does bias not go up when covariates omitted?
## I guess because these covariates are not strongly associated with transition rates
## ollwsdrt is more associated, particularly with 2-3 rate, so bigger biases when ollwsdrt effect omitted 
## Variance seems to mostly go down as covariates omitted 

pp <- 3 # intercepts always included
fic.msm(psor.wide.msm, inds=c(0,1,1,1,1,1), pp, focus)
fic.msm(psor.wide.msm, inds=c(0,0,1,1,1,1), pp, focus)
fic.msm(psor.wide.msm, inds=c(0,0,0,1,1,1), pp, focus)
fic.msm(psor.wide.msm, inds=c(0,0,0,0,1,1), pp, focus)
fic.msm(psor.wide.msm, inds=c(0,0,0,0,0,1), pp, focus)


}


## function below is included from 2017-12-08 in msm development version on
## https://cran.r-project.org/package=msm

## updatepars.msm <- function(x, pars){
##     x.rep <- x
##     x.rep$paramdata$params <- pars
##     x.rep$estimates <- pars
##     output <- msm.form.output(x.rep, "intens")
##     x.rep$Qmatrices <- output$Qmatrices
##     if (x$emodel$misc) {
##         output <- msm.form.output(x.rep, "misc")
##         x.rep$Ematrices <- output$Ematrices
##         names(x.rep$Ematrices)[1] <- "logitbaseline"
##     }
##     x.rep    
## }
