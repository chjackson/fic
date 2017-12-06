# library(devtools)
# load_all("..")
# document("..")
# library(MASS) 

birthwt <- within(birthwt, { 
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
X = with(birthwt, cbind(intercpt,lwt))
Z = with(birthwt, cbind(age, smoke, ht,ui, smokeage=age*smoke, smokeui=smoke*ui))


###### Compute FIC using new code

wide.glm <- glm(Y ~ X + Z - 1, data=birthwt, family=binomial)
ests <- coef(wide.glm)
n <- nobs(wide.glm)
J <- solve(vcov(wide.glm)) / n
pp <- ncol(X)
## Focus quantity: 
## Prob of LBW for covariate values of the first person in the data
vals.first <- c(X[1,], Z[1,])
focus <- function(ests){
    plogis(q = ests %*% vals.first)
}
## Submodel to compare with the wide model: all covariates included except for the last 
inds <- c(1,1,1,1,1,0)

fic(ests, J, inds, pp, n, focus)


###### Compare against FIC using Gerda's original code

source("FIClogisticreg.R")
ficall <- FIC.logistic.regression(Y=Y, X=X, Z=Z, XZeval=cbind(X,Z), dataframe=NULL)
## extract results for focus: LBW for covariate values of first person in data
## gives matrix with one row for each submodel 
res <- sapply(ficall[c("FIC","Bias","Bias2","Var","VarS")], function(x)x[,1])
## row 2 is submodel 1,1,1,1,1,0
res[2,]  ## matches my code to small d.p, differences due to numerical derivatives? 
