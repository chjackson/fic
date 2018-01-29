## ------------------------------------------------------------------------
library(msm)
psor.q <- rbind(c(0,0.1,0,0), c(0,0,0.1,0), c(0,0,0,0.1), c(0,0,0,0))
psor.wide.msm <- msm::msm(state ~ months, subject=ptnum, data=msm::psor, qmatrix = psor.q,  covariates = ~ollwsdrt+hieffusn, control=list(fnscale=1))
psor.wide.msm

## ------------------------------------------------------------------------
#install_github("chjackson/msm") # if necessary 
focus_tlos <- function(pars){
    x.new <- msm::updatepars.msm(psor.wide.msm, pars)
    msm::totlos.msm(x.new, covariates=0, t=10)["State 4"]
}

## ------------------------------------------------------------------------
inds0 <- c(1,1,1,0,0,0,0,0,0) # intercepts always included
inds <- rbind(
    c(1,1,1,1,1,1,1,1,1),
    c(1,1,1,0,1,1,1,1,1),
    c(1,1,1,0,0,1,1,1,1),
    c(1,1,1,0,0,0,1,1,1),
    c(1,1,1,0,0,0,0,1,1),
    c(1,1,1,0,0,0,0,0,1),
    c(1,1,1,0,0,0,0,0,0))
fic(wide=psor.wide.msm, inds=inds, inds0=inds0, focus=focus_tlos)

