---
title: "Using the fic R package for focused model comparison: multi-state models"
date: "2018-06-18"
author: Chris Jackson <chris.jackson@mrc-bsu.cam.ac.uk>
output: 
  html_document: 
    includes:
      before_body: fic.sty  
  code_folding: show
vignette: >
 %\VignetteIndexEntry{Using the fic R package for focused model comparison: multi-state models}
 %\VignetteEngine{knitr::rmarkdown}
 %\VignetteDepends{fic,ggplot2,msm}
 %\VignetteEncoding{UTF-8}
bibliography: fic.bib
---

## Focused model comparison for covariate selection in multi-state models fitted with msm 

Example dataset in package, `psor`.  Four state progression-only model.

Wide model: model with two binary covariates associated with different effects for all three transition rates.


```r
if (!require("msm")) stop("The `msm` package should be installed to run code in this vignette") 
```

```
## Loading required package: msm
```

```r
psor.q <- rbind(c(0,0.1,0,0), c(0,0,0.1,0), c(0,0,0,0.1), c(0,0,0,0))
psor.wide.msm <- msm(state ~ months, subject=ptnum, data=psor, qmatrix = psor.q,  covariates = ~ollwsdrt+hieffusn, control=list(fnscale=1))
psor.wide.msm
```

```
## 
## Call:
## msm(formula = state ~ months, subject = ptnum, data = psor, qmatrix = psor.q,     covariates = ~ollwsdrt + hieffusn, control = list(fnscale = 1))
## 
## Maximum likelihood estimates
## Baselines are with covariates set to their means
## 
## Transition intensities with hazard ratios for each covariate
##                   Baseline                    ollwsdrt              
## State 1 - State 1 -0.1004 (-0.12750,-0.07898)                       
## State 1 - State 2  0.1004 ( 0.07898, 0.12750) 0.7320 (0.4258,1.2585)
## State 2 - State 2 -0.1623 (-0.20601,-0.12789)                       
## State 2 - State 3  0.1623 ( 0.12789, 0.20601) 0.4579 (0.2643,0.7932)
## State 3 - State 3 -0.2607 (-0.34952,-0.19453)                       
## State 3 - State 4  0.2607 ( 0.19453, 0.34952) 1.5757 (0.7776,3.1928)
##                   hieffusn            
## State 1 - State 1                     
## State 1 - State 2 2.338 (1.0937,4.997)
## State 2 - State 2                     
## State 2 - State 3 1.681 (0.9500,2.975)
## State 3 - State 3                     
## State 3 - State 4 1.394 (0.7738,2.511)
## 
## -2 * log-likelihood:  1112.613 
## [Note, to obtain old print format, use "printold.msm"]
```

This requires version 1.6.6 of `msm`, available from CRAN since 3 Feb 2017.  This version introduced the `updatepars.msm` function for altering the point estimates from a fitted model a model to a vector of values supplied by the user.  This allows functions such as `totlos.msm` to be used, which define complicated functions of the model parameters.

Here we use this to define a focus function, which returns the expected total time spent in state 4 over 10 years for people without `ollwsdrt` or `hieffusn`, given a vector of parameters `pars` in the model structure `psor.wide.msm`. 


```r
focus_tlos <- function(pars){
    x.new <- updatepars.msm(psor.wide.msm, pars)
    totlos.msm(x.new, covariates=0, t=10)["State 4"]
}
```

We assess the wide model and six further submodels.



