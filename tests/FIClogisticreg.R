### Functions from
### https://perswww.kuleuven.be/~u0043181/modelselection/index.html#Programs_for_the_focussed_information


combinations = function(n){
  comb = NULL
  if (n<15) {
  for( i in 1:n) comb = rbind(cbind(1,comb),cbind(0,comb))
  return(comb)
  }
  else {error("this value will probably hang your computer, try on your own risk")}
 }

# Example, models with extra variables 2 and 3
#variables = c(0,1,1,0,0)

FIC.main = function(variables,K,omega,deltahat,psi.full)
  {
    if (sum(variables)>0)
      {
        qq = length(variables)
        indic = variables*(1:qq)
        Id =  diag(rep(1,qq))
        Kinv = solve(K)

        pi.S = matrix(Id[indic,],ncol=qq)

        K.S = solve(pi.S %*% Kinv %*% t(pi.S))
	  M.S = t(pi.S) %*% K.S %*% pi.S %*% Kinv

        psi.S = t(omega)%*% M.S %*% deltahat

        omega.S = pi.S %*% omega

        FIC.S = (psi.full-psi.S)^2 +2*t(omega.S)%*%K.S%*%omega.S

        bias.S = t(omega)%*%deltahat-psi.S

        bias2.S = sign(bias.S)*sqrt(max(0,t(omega)%*%(Id-M.S)%*%(deltahat%*%t(deltahat)-K)%*%(Id-M.S)%*%omega))
    
	  var2.S = t(omega) %*% M.S %*% omega 
	  }
    else {FIC.S = psi.full^2
    bias.S = bias2.S = t(omega) %*% deltahat
    var2.S = 0}
    
    outputFIC=list(FIC.S,bias.S,bias2.S,var2.S)
    names(outputFIC)=c("FIC.S","bias.S","bias2.S","var2.S")
    return(outputFIC) 
  }


invlogit =  function(u){exp(u) / (1 + exp(u))} 

FIC.logistic.regression = function(Y,X,Z,XZeval=cbind(X,Z),dataframe)
{
  qq = ncol(Z)
  pp = ncol(X) 
  n = nrow(X)
  XZ=cbind(X,Z)  

  # Fit the full model
  fit.full =  glm(Y ~ X + Z -1, family=binomial)  # intercept already in X matrix
  betagamma.hat = fit.full$coef 
  est.prob = fit.full$fitted.values

  J <- matrix(data=NA, nrow=(pp+qq), ncol=(pp+qq))
  for (i in 1:(pp+qq))
   {for (j in 1:(pp+qq))
      {J[i,j] <- mean(est.prob*(1-est.prob)*XZ[ ,i]*XZ[ ,j])}}

  deltahat = n^{1/2}*betagamma.hat[-(1:pp)]
  J00 <- J[1:pp, 1:pp]
  J10 <- J[-(1:pp), 1:pp]
  J01 <- J[1:pp,-(1:pp)]
  J11 <- J[-(1:pp),-(1:pp)]

  invJ = solve(J)
  K <- invJ[-(1:pp),-(1:pp)]


  # Make all combinations of 0/1 to form all subsets of the full model
  varmatrix = combinations(qq)


  ## Main part of the program!
  ############################

  # This matrix will contain per subject the FIC value for each of the models.

  n.eval = nrow(XZeval)
  risk.model.subject = FIC.model.subject = matrix(NA,ncol=n.eval,nrow=nrow(varmatrix))
  bias.model.subject = bias2.model.subject = var.model.subject = risk.model.subject
  varS.model.subject = muS = var.model.subject
  aic.model = bic.model = rep(NA,nrow(varmatrix))

  for (k in 1:n.eval)
    {
      # print(k)
      x.eval = as.numeric(XZeval[k,1:pp])
      z.eval = as.numeric(XZeval[k,-(1:pp)])
      
      p0 = invlogit(t(betagamma.hat)%*%c(x.eval,z.eval))
      omega = (J10%*%solve(J00)%*%x.eval-z.eval)*as.vector(p0*(1-p0))

      psi.full = t(omega) %*% deltahat

	partmutheta = x.eval*as.vector(p0*(1-p0))
	tau0sq = t(partmutheta)%*%solve(J00)%*%partmutheta

      for (j in 1:nrow(varmatrix)) 
	{
         FIC.out = FIC.main(varmatrix[j,],K=K,omega=omega,
         deltahat=deltahat,psi.full=psi.full)

         FIC.model.subject[j,k]=FIC.out$FIC.S
	   #write(FIC.model.subject[,k],file=outfile,ncol=nrow(varmatrix),append=T)

         risk.model.subject[j,k] = tau0sq-t(omega)%*%K%*%omega+FIC.model.subject[j,k]
         bias.model.subject[j,k] = FIC.out$bias.S/n
         bias2.model.subject[j,k] = FIC.out$bias2.S/n
	   var.model.subject[j,k] = (tau0sq + FIC.out$var2.S)/n
         calc.in.S = Inside.model.S(varmatrix[j,],Y,X,Z,x.eval,z.eval)
         varS.model.subject[j,k] = calc.in.S$var/n
         muS[j,k] = calc.in.S$mu.est
         if (k==1) { aic.model[j] = calc.in.S$aic; bic.model[j] = calc.in.S$bic }
       } 
    NULL
    } # end loop over k

  risksqrt = sqrt(risk.model.subject/n)  
  
  return(list(Est=muS,FIC=FIC.model.subject,Risksqrt=risksqrt,
		Bias=bias.model.subject,Bias2=bias2.model.subject,Var=var.model.subject,
		VarS=varS.model.subject,AIC=aic.model,BIC=bic.model))
}


Inside.model.S = function(variables,YY,XX,ZZ,xeval,zeval)
{
  qq = length(variables)
  pp = ncol(XX)
  indic = variables*(1:qq)

  if (sum(variables)>0)
  {     XZ.S = cbind(XX,ZZ[,indic])
        fit.S =  glm(YY ~ XX + ZZ[,indic] -1, family=binomial)  
	  xzeval = c(xeval,zeval[indic])
  }
  else 
  {     XZ.S = XX
        fit.S =  glm(YY ~ XX -1, family=binomial)  
	  xzeval = xeval
  } 
  betagamma.hat.S = fit.S$coef 
  est.prob.S = fit.S$fitted.values

  # J.S : fill in estimators of model S itself
  pq.S = ncol(XZ.S)
  J.S = matrix(data=NA, nrow=(pq.S), ncol=(pq.S))
  for (i in 1:(pq.S))
  	{for (j in 1:(pq.S))
      	{  J.S[i,j] = mean(est.prob.S*(1-est.prob.S)*XZ.S[ ,i]*XZ.S[ ,j]) }}

  p0.S = invlogit(t(betagamma.hat.S)%*%xzeval)
  J00.S <- J.S[1:pp, 1:pp]
  if (sum(variables)>0)
  {     J10.S <- J.S[-(1:pp), 1:pp]
        J01.S <- J.S[1:pp,-(1:pp)]
        J11.S <- J.S[-(1:pp),-(1:pp)]
        omega.S = (J10.S%*%solve(J00.S)%*%xeval-zeval[indic])*as.vector(p0.S*(1-p0.S))
  }
  else  omega.S=0
  partmutheta.S = xeval*as.vector(p0.S*(1-p0.S))
  tau0sq.S = t(partmutheta.S)%*%solve(J00.S)%*%partmutheta.S
  var.model.subject.S = (tau0sq.S + t(omega.S) %*% omega.S)

  bic.S = -AIC(fit.S,k=log(length(Y)))
  aic.S = -AIC(fit.S,k=2)
  
  return(list(mu.est=p0.S,aic=aic.S,bic=bic.S,var=var.model.subject.S))
} 
