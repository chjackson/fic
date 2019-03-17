## ------------------------------------------------------------------------
if (!require("msm"))
    stop("The `msm` package should be installed
to run code in this vignette") 
Qind <- rbind(c(0, 1, 0, 0),
              c(0, 0, 1, 0),
              c(0, 0, 0, 1),
              c(0, 0, 0, 0))
psor.wide.msm <- msm(state ~ months, subject=ptnum, data=psor, 
                     qmatrix = Qind,  gen.inits=TRUE,
                     covariates = ~ollwsdrt+hieffusn)
psor.wide.msm

## ------------------------------------------------------------------------
inds <- rbind(
    c(1,1,1,0,0,0,0,0,0),
    c(1,1,1,0,0,0,0,0,1),
    c(1,1,1,0,0,0,0,1,1),
    c(1,1,1,0,0,0,1,1,1),
    c(1,1,1,0,0,1,1,1,1),
    c(1,1,1,0,1,1,1,1,1),
    c(1,1,1,1,1,1,1,1,1)
)

## ------------------------------------------------------------------------
totlos.msm(psor.wide.msm, covariates=0, tot=10)

## ------------------------------------------------------------------------
focus_tlos <- function(par){
    x.new <- updatepars.msm(psor.wide.msm, par)
    totlos.msm(x.new, covariates=0, tot=10)["State 4"]
}

## ------------------------------------------------------------------------
library(fic)
fic(wide=psor.wide.msm, inds=inds, focus=focus_tlos)

