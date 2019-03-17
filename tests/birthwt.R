# library(devtools)
# library(testthat)
# load_all("..")
# document("..")

# processed birthwt data included in fic package

# full Z matrix 
#Z = with(birthwt, cbind(age,raceblack,raceother,smoke,ptl,ht,ui,ftv1,ftv2p,
#                        smokeui=smoke*ui,ageftv1=age*ftv1,ageftv2=age*ftv2p,smokeage=age*smoke))

## Example 6.1 in book 
## For the present illustration, we include in every candidate model the intercept x 1 = 1 and the weight x 2 (in kg) of the mother prior to pregnancy. Other covariates, from which we wish to select a relevant subset, are [ age, smoking, hypertension, uterine irritability, smok*age, smok*ut ] 

Y = birthwt$low
X = with(birthwt, cbind(intercpt, lwtkg))
Z = with(birthwt, cbind(age, smoke, ht, ui, smokeage, smokeui))
XZ <- cbind(X, Z)
    
###### Compute FIC using new code

wide.glm <- glm(Y ~ XZ - 1, data=birthwt, family=binomial)
## Submodel to compare with the wide model: all covariates included except for the last 
inds <- c(1,1,1,1,1,1,1,0)
sub.glm <- glm(Y ~ XZ[,inds==1] - 1, data=birthwt, family=binomial)
## First two covariates always included
inds0 <- c(1,1,0,0,0,0,0,0)
## Focus quantity: 
## Prob of LBW for covariate values of the first person in the data
vals.first <- XZ[1,]
focus <- function(par){
    plogis(q = par %*% vals.first)
}

fic(wide=wide.glm, inds=inds, inds0=inds0, focus=focus)

pred <- predict.glm(wide.glm, type="response", se.fit=TRUE)
cbind(pred$fit, pred$se.fit)[1,]
fic(wide=wide.glm, inds=c(1,1,1,1,1,1,1,1), inds0=inds0, focus=focus) 
## check SE of focus under wide model matches RMSE returned by FIC.  
## OK - this is tau0sq in code, which matches the delta method formula


###### Compare against FIC using Gerda's original code

source("FIClogisticreg.R")
ficall <- FIC.logistic.regression(Y=Y, X=X, Z=Z, XZeval=XZ[1:2,], dataframe=NULL)
## extract results for focus: LBW for covariate values of first and second people in data
## gives matrix with one row for each submodel 
do.call("rbind", ficall[c("FIC","Bias","Bias2","Var","VarS")])
## matches my code to first few sig. fig. 


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
fic.ana <- fic(wide=wide.glm, inds=inds, inds0=inds0, focus="prob_logistic", X=vals.first)
fic.ana

fic.num <- fic_multi(par=coef(wide.glm), J=solve(vcov(wide.glm))/nobs(wide.glm), inds=inds, inds0=inds0, n=nobs(wide.glm), focus=focus)
fic.num

library(testthat)
expect_equal(fic.ana[,"FIC"], fic.num[,"FIC",1]) # matches within reasonable precision 

## Alternatively: supply derivative values directly to lower-level fic() function 
object <- wide.glm
ests <- coef(object)
n <- nobs(object)
J <- solve(vcov(object)) / n
focus_deriv <- t(prob_logistic_deriv(ests, vals.first))
fic.ana2 <- fic_core(par=ests, J=J, inds=inds, inds0=inds0, n=n, focus_deriv=focus_deriv, parsub=coef(sub.glm))
expect_equal(as.numeric(fic.ana2[,"FIC"]), fic.ana[,"FIC"])

