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
inds <- rbind(mod1 = c(1,1,1,1,0,0,0,0),
              mod2 = c(1,1,1,1,1,0,0,0))
inds0 <- c(1,1,0,0,0,0,0,0)
sub <- list(mod1.glm, mod2.glm)
fic1 <- fic(wide=wide.glm, inds=inds, inds0=inds0, focus=focus, X=X, sub=sub)
fic1

## ------------------------------------------------------------------------
fns <- list(coef = function(x)coef(x),
            nobs = function(x)nobs(x),
            vcov = function(x)vcov(x))
fic1 <- fic(wide=wide.glm, inds=inds, inds0=inds0, focus=focus, fns=fns, X=X, sub=sub)

## ------------------------------------------------------------------------
combs <- expand.grid(intercept=1, lwtkg=1, age=c(0,1), smoke=c(0,1), ht=c(0,1), ui=c(0,1), smokeage=c(0,1), smokeui=c(0,1))

## ------------------------------------------------------------------------
ficres <- fic(wide=wide.glm, inds=combs, inds0=inds0, focus=focus, X=X)

## ------------------------------------------------------------------------
combs <- with(combs,
              combs[!((smoke==0 & smokeage==1) |
                      (smoke==0 & smokeui==1) |
                      (age==0 & smokeage==1) |
                      (ui==0 & smokeui==1)),])
rownames(combs) <- apply(combs, 1, paste, collapse="")
namesfull <- apply(combs, 1, function(x)paste(colnames(combs)[x==1], collapse=","))
ficres <- fic(wide=wide.glm, inds=combs, inds0=inds0, focus=focus, X=X)

## ----warning=FALSE-------------------------------------------------------
Xobs <- with(birthwt, cbind(intercpt, lwtkg))
Zobs <- with(birthwt, cbind(age, smoke, ht, ui, smokeage, smokeui))
nmod <- nrow(combs)
sub <- vector(nmod, mode="list")
for (i in 1:nmod){
  XZi <- cbind(Xobs, Zobs)[,which(combs[i,]==1)]
  sub[[i]] <- glm(low ~ XZi - 1, data=birthwt, family=binomial)
}
ficres <- fic(wide=wide.glm, inds=combs, inds0=inds0, focus=focus, X=X, sub=sub)
ficres

## ----message=FALSE-------------------------------------------------------
library(ggplot2)

## ------------------------------------------------------------------------
ress <- within(ficres,{
  l95 = focus - qnorm(0.975)*se
  u95 = focus + qnorm(0.975)*se
})
ress <- subset(ress, !is.nan(ress$rmse))
ress <- ress[order(ress$rmse, decreasing=TRUE),]
ress_wide <- ress[ress$mods=="11111111",]
wide_est <- focus(coef(wide.glm), X)
ress$wide_est <- wide_est[match(ress$vals, rownames(wide_est))]

ps <- ggplot(ress[ress$mods!="11111111",],
         aes(x=focus, y=rmse)) +
  facet_grid(.~vals) + 
  xlim(0, 0.6) + 
  xlab("Probability of low birth weight") +
  ylab("Root mean square error of estimate") + 
  geom_point() +
  geom_segment(aes(x=l95, xend=u95, yend=rmse)) +
  geom_point(data=ress_wide, col="red") +
  geom_segment(aes(x=l95, xend=u95, yend=rmse), data=ress_wide, col="red") +
  geom_text(aes(x=0, label="11111111", hjust=0), data=ress_wide, col="red", size=5) +
  geom_text(aes(x=0, label=mods, hjust=0)) +
  geom_vline(aes(xintercept = wide_est), col="red") +
  theme(text = element_text(size=20),
        axis.text.x = element_text(size=20)) +
  geom_text(aes(x=wide_est, y= 0.12, label="Estimate from wide model", hjust=0), col="red") 
ps

