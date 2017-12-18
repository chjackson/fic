# library(devtools)
# library(testthat)
# load_all("..")
# document("..")
# library(MASS) 

if (requireNamespace("MASS", quietly = TRUE)){

birthwt <- within(MASS::birthwt, { 
    lwtkg = lwt*0.45359237 # weight measured in lb in MASS data, kg in FIC book
    race = factor(race, labels = c("white", "black", "other"))
    ptd = factor(ptl > 0)
    ftv = factor(ftv) # keeping all levels, not transforming to 2+
    levels(ftv)[-(1:2)] <- "2+" # doing the transformation
    ftv1=as.numeric(ftv=="1")
    ftv2p=as.numeric(ftv=="2+")
    low = factor(low)
    smoke = (smoke > 0)
    ht = (ht > 0)
    ui = (ui > 0)
    raceblack <- 1*(race=="black")
    raceother <- 1*(race=="other")
    intercpt = 1 + 0*seq(1:length(low)) # intercept
}
)
# full Z matrix 
#Z = with(birthwt, cbind(age,raceblack,raceother,smoke,ptl,ht,ui,ftv1,ftv2p,
#                        smokeui=smoke*ui,ageftv1=age*ftv1,ageftv2=age*ftv2p,smokeage=age*smoke))

## Example 6.1 in book 
## For the present illustration, we include in every candidate model the intercept x 1 = 1 and the weight x 2 (in kg) of the mother prior to pregnancy. Other covariates, from which we wish to select a relevant subset, are [ age, smoking, hypertension, uterine irritability, smok*age, smok*ut ] 

Y = birthwt$low
X = with(birthwt, cbind(intercpt,lwtkg))
Z = with(birthwt, cbind(age, smoke, ht,ui, smokeage=age*smoke, smokeui=smoke*ui))
XZ <- cbind(X, Z)
    
###### Compute FIC using new code

wide.glm <- glm(Y ~ X + Z - 1, data=birthwt, family=binomial)
## Submodel to compare with the wide model: all covariates included except for the last 
inds <- c(1,1,1,1,1,0)
## First two covariates always included
pp <- ncol(X)
## Focus quantity: 
## Prob of LBW for covariate values of the first person in the data
vals.first <- c(X[1,], Z[1,])
focus <- function(ests){
    plogis(q = ests %*% vals.first)
}

fic.glm(wide.glm, inds, pp, focus)

pred <- predict.glm(wide.glm, type="response", se.fit=TRUE)
cbind(pred$fit, pred$se.fit)[1,]
fic.glm(wide.glm, inds=c(1,1,1,1,1,1), pp, focus) # check SE of focus under wide model matches RMSE returned by FIC.  OK - this is tau0sq in code, which matches the delta method formula


###### Compare against FIC using Gerda's original code

source("FIClogisticreg.R")
ficall <- FIC.logistic.regression(Y=Y, X=X, Z=Z, XZeval=XZ[1:2,], dataframe=NULL)
## extract results for focus: LBW for covariate values of first person in data
## gives matrix with one row for each submodel 
res <- sapply(ficall[c("FIC","Bias","Bias2","Var","VarS")], function(x)x[,1])
## row 2 is submodel 1,1,1,1,1,0
res[2,]  ## matches my code to first few sig. fig. 


    ## Differences arise from differences in the information matrix J.
    ## My code uses vcov() on the fitted GLM object.
    ## Gerda's code calculates observed information matrix by hand, evaluating second derivs of the log likelihood at the MLEs
    ## so where does the estimate used by vcov.glm() come from?
    # vcov calculated as follows: 
    # W <- wide.glm$weights
    # solve(t(XZ) %*% (W * XZ)) # yes. this matches vcov
    #p0 <- wide.glm$fitted.values
    #cbind(p0*(1-p0), wide.glm$weights) # small differences here.
    # Venables and Ripley p186 implies that this is the expected, not the observed information. presumably more details in Mccullagh and Nelder.

    
    ### Compare numerical and analytic derivatives
    fic.ana <- fic.glm(wide.glm, inds, pp, focus="prob_logistic", X=vals.first)
    fic.ana
    
    fic.num <- fic(ests=ests, J=J, inds=inds, pp=pp, n=n, focus=focus)
    expect_equal(fic.ana, fic.num) # matches within reasonable precision 

    ## Alternatively: supply derivative values directly to lower-level fic() function 
    object <- wide.glm
    ests <- coef(object)
    n <- nobs(object)
    J <- solve(vcov(object)) / n
    focus_deriv <- prob_logistic_deriv(ests, vals.first)
    fic.ana2 <- fic(ests=ests, J=J, inds=inds, pp=pp, n=n, focus_deriv=focus_deriv)
    expect_equal(fic.ana2, fic.ana)
    
    ## test error handling 
    expect_error(fic.glm(wide.glm, inds, pp, focus="nonexistent_function", X=vals.first),
                 "not found")

}
