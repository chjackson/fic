#library(devtools)
#load_all("../../../../msm/msm/")
#library(msm)
#load_all("..")
#install_github("chjackson/msm") # requires devel version 

if (requireNamespace("msm", quietly = TRUE)){

## 4-state PsA model using "psor" data in msm
## Wide model includes two covariates on all rates 
psor.q <- rbind(c(0,0.1,0,0), c(0,0,0.1,0), c(0,0,0,0.1), c(0,0,0,0))
psor.wide.msm <- msm::msm(state ~ months, subject=ptnum, data=msm::psor, qmatrix = psor.q,  covariates = ~ollwsdrt+hieffusn, control=list(fnscale=1))
psor.wide.msm

## focus: total length of stay in state 4 up to 10 years for people without ollwsdrt or hieffusn
focus_tlos <- function(pars){
    x.new <- msm::updatepars.msm(psor.wide.msm, pars)
    msm::totlos.msm(x.new, covariates=0, t=10)["State 4"]
}

msm::totlos.msm(psor.wide.msm, covariates=0, t=10, ci="normal") # SE around 0.1
focus_tlos(psor.wide.msm$estimates) # should match 

## FIC analysis: omit covariate effects in turn from full model 
## Why does bias not always go up when covariates omitted?
## I guess because these covariates are not strongly associated with transition rates
## ollwsdrt is more associated, particularly with 2-3 rate, so (broadly) bigger biases when ollwsdrt effect omitted, though seems to be noise in the bias estimates. 
## Variance does go down smoothly as covariates omitted

inds0 <- c(1,1,1,0,0,0,0,0,0) # intercepts always included
fic.msm(wide=psor.wide.msm, inds=c(1,1,1,1,1,1,1,1,1), inds0=inds0, focus=focus_tlos)
fic.msm(wide=psor.wide.msm, inds=c(1,1,1,0,1,1,1,1,1), inds0=inds0, focus=focus_tlos)
fic.msm(wide=psor.wide.msm, inds=c(1,1,1,0,0,1,1,1,1), inds0=inds0, focus=focus_tlos)
fic.msm(wide=psor.wide.msm, inds=c(1,1,1,0,0,0,1,1,1), inds0=inds0, focus=focus_tlos)
fic.msm(wide=psor.wide.msm, inds=c(1,1,1,0,0,0,0,1,1), inds0=inds0, focus=focus_tlos)
fic.msm(wide=psor.wide.msm, inds=c(1,1,1,0,0,0,0,0,1), inds0=inds0, focus=focus_tlos)
fic.msm(wide=psor.wide.msm, inds=c(1,1,1,0,0,0,0,0,0), inds0=inds0, focus=focus_tlos)


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
