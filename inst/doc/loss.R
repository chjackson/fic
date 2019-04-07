## ------------------------------------------------------------------------
library(fic)
wide.glm <- glm(low ~ lwtkg + age + smoke + ht + ui + smokeage + smokeui, 
                data=birthwt, family=binomial)
vals.smoke <- c(1, 58.24, 22.95, 1, 0, 0, 22.95, 0)
X <- vals.smoke
inds0 <- c(1,1,0,0,0,0,0,0)
combs <- all_inds(wide.glm, inds0)
ficres <- fic(wide=wide.glm, inds=combs, inds0=inds0, 
              focus=prob_logistic, X=X)

## ------------------------------------------------------------------------
set.seed(1)
ficboot_mse <- fic(wide=wide.glm, inds=combs, inds0=inds0, 
                   focus=prob_logistic, X=X, B=1000)

## ----fig.width=5,fig.height=5--------------------------------------------
plot(ficres$rmse, ficboot_mse$loss, xlim=c(0.05,0.15), 
     ylim=c(0.05,0.15), pch=19, 
     xlab = "Root mean square error under FIC asymptotic theory",
     ylab = "Root mean square error from bootstrap under wide model")
abline(a=0, b=1)

## ------------------------------------------------------------------------
loss_abserror <- function(sub, wide){
    mean(abs(sub - wide))
}
ficboot_abs <- fic(wide=wide.glm, inds=combs, inds0=inds0, 
                   focus=prob_logistic, X=X, B=1000, loss=loss_abserror)

