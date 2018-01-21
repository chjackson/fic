## ------------------------------------------------------------------------
library(fic)
wide.glm <- glm(low ~ lwtkg + age + smoke + ht + ui + smokeage + smokeui, data=birthwt, family=binomial)

## ------------------------------------------------------------------------
focus <- function(par, X)plogis(X %*% par)
vals.smoke <-    c(1, 58.24, 22.95, 1, 0, 0, 22.95, 0)
vals.nonsmoke <- c(1, 59.50, 23.43, 0, 0, 0, 0, 0)
X <- rbind(vals.smoke, vals.nonsmoke)

## ------------------------------------------------------------------------
focus(coef(wide.glm), X=X)

## ------------------------------------------------------------------------
mod1.glm <- glm(low ~ lwtkg + age + smoke, data=birthwt, family=binomial)
mod2.glm <- glm(low ~ lwtkg + age + smoke + ht, data=birthwt, family=binomial)

## ------------------------------------------------------------------------
inds0 <- c(1,1,0,0,0,0,0,0)
inds1 <- c(1,1,1,1,0,0,0,0)
inds2 <- c(1,1,1,1,1,0,0,0)
fic1 <- fic(wide=wide.glm, sub=mod1.glm, inds=inds1, inds0=inds0, focus=focus, X=X)
fic2 <- fic(wide=wide.glm, sub=mod2.glm, inds=inds2, inds0=inds0, focus=focus, X=X)
fic1
fic2

## ------------------------------------------------------------------------
par <- coef(wide.glm)
J <- solve(vcov(wide.glm))/nrow(birthwt)
fic(par=par, J=J,  inds=inds1, inds0=inds0, n=nrow(birthwt), focus=focus, X=X, parsub=coef(mod1.glm))

## ------------------------------------------------------------------------
combs <- expand.grid(age=c(0,1), smoke=c(0,1), ht=c(0,1), ui=c(0,1), smokeage=c(0,1), smokeui=c(0,1))
combs <- as.matrix(with(combs,
                        combs[!((smoke==0 & smokeage==1) |
                                (smoke==0 & smokeui==1) |
                                (age==0 & smokeage==1) |
                                (ui==0 & smokeui==1)),]))
rownames(combs) <- apply(combs, 1, paste, collapse="")
namesfull <- apply(combs, 1, function(x)paste(colnames(combs)[x==1], collapse=","))

## ----warning=FALSE-------------------------------------------------------
Xobs <- with(birthwt, cbind(intercpt, lwtkg))
Zobs <- with(birthwt, cbind(age, smoke, ht, ui, smokeage, smokeui))

res <- array(dim=c(nrow(X), nrow(combs), 9))
for (i in 1:nrow(combs)){
  XZi <- cbind(Xobs, Zobs[,which(combs[i,]==1)])
  sub.glm <- glm(low ~ XZi - 1, data=birthwt, family=binomial)
  ficres <- fic(wide=wide.glm, sub=sub.glm, 
                inds=c(1,1,combs[i,]), inds0=inds0, focus=focus, X=X)
  res[,i,] <- ficres
}
dimnames(res) <- list(rownames(X), rownames(combs), colnames(ficres))

## ----message=FALSE-------------------------------------------------------
library(tidyverse)

## ------------------------------------------------------------------------
ress <- as.data.frame(res["vals.smoke",,]) %>% 
  mutate(l95 = focus - qnorm(0.975)*se,
         u95 = focus + qnorm(0.975)*se,
         namesfull = namesfull) %>% 
  filter(!is.nan(rmse)) %>%
  arrange(desc(rmse)) %>%
  mutate(aic.rank = rank(AIC),
         bic.rank = rank(BIC))

mod_wide <-  namesfull[rowSums(combs)==ncol(combs)]
ress_wide <- ress[ress$namesfull==mod_wide,]
mods_exc_smoke <- c(mod_wide, 
                    "smoke",
                    "smoke,ui,smokeui",
                    "age,smoke,ht,ui,smokeui",
                    "smoke,ht,ui",
                    "smoke,ui",
                    "smoke,ht",
                    "age,smoke,ht"
                    )

ps <- ress %>% filter(!namesfull %in% mods_exc_smoke) %>% 
  ggplot(aes(x=focus, y=rmse)) +
  xlim(0, 0.6) + 
  xlab("Probability of low birth weight") +
  ylab("Root mean square error of estimate") + 
  geom_point() +
  geom_segment(aes(x=l95, xend=u95, yend=rmse)) +
  geom_point(data=ress_wide, col="red") +
  geom_segment(aes(x=l95, xend=u95, yend=rmse), data=ress_wide, col="red") +
  geom_text(aes(x=0, label=mod_wide, hjust=0), data=ress_wide, col="red", size=5) +
  #  geom_point(aes(x = focus - bias), col="red") +
  geom_text(aes(x=0, label=namesfull, hjust=0)) +
  geom_vline(aes(xintercept = as.numeric(plogis(coef(wide.glm) %*% vals.smoke))),
             col="red") +
  theme(text = element_text(size=20),
        axis.text.x = element_text(size=20)) +
  geom_text(aes(x=0.35, y= 0.08, label="Estimate from wide model", hjust=0), col="red") +
  geom_text(aes(x=0, y=0.047, label="Models with lowest MSE", hjust=0), col="blue") +
  geom_text(aes(x=0, y=0.085, label="Models with highest MSE", hjust=0), col="blue") +
  geom_text(aes(x=0.55, label=aic.rank), size=4) + 
  geom_text(data=ress_wide, aes(x=0.55, label=aic.rank), size=5, col="red") +
  geom_text(aes(x=0.57, y=0.093, label="Ranks of:")) + 
  geom_text(aes(x=0.55, y=0.091, label="AIC")) +
  geom_text(aes(x=0.6, label=bic.rank), size=4) + 
  geom_text(data=ress_wide, aes(x=0.6, label=bic.rank), size=5, col="red") +
  geom_text(aes(x=0.6, y=0.091, label="BIC")) 
ps

