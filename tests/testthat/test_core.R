par <- coef(wide.glm)
J <- solve(vcov(wide.glm))/nrow(birthwt)
fic_core(par=par, J=J,  inds=inds1, inds0=inds0, n=nrow(birthwt), focus=focus, X=X, parsub=coef(mod1.glm))
