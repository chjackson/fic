## ----echo=FALSE-----------------------------------------------------
options(prompt = "R> ",
        continue = "+  ",
        width = 70,
        useFancyQuotes = FALSE,
        digits = 3)
library("knitr")
opts_chunk$set(fig.path="fic-")

## ----warning=FALSE,message=FALSE------------------------------------
library(GGally)
mtcars$am <- factor(mtcars$am)
ggpairs(mtcars[,c("mpg","am","wt","qsec","disp","hp")], aes(colour=am))

## -------------------------------------------------------------------
wide.lm <- lm(mpg ~ am + wt + qsec + disp + hp, data=mtcars)

## -------------------------------------------------------------------
library(fic)
ncovs_wide <- length(coef(wide.lm)) - 1
inds0 <- c(1, rep(0, ncovs_wide))
inds <- all_inds(wide.lm, inds0)

## -------------------------------------------------------------------
cmeans <- colMeans(model.frame(wide.lm)[,c("wt","qsec","disp","hp")])
X <- rbind(
  "auto"   = c(intercept=1, am=0, cmeans),
  "manual" = c(intercept=1, am=1, cmeans)
)
ficres <- fic(wide.lm, inds=inds, focus=mean_normal, X=X)
summary(ficres)

## ----mtcars,include=FALSE-------------------------------------------
ggplot_fic(ficres)

## -------------------------------------------------------------------
library(gapminder)
gap2 <- gapminder[gapminder$country !="Kuwait",]
pal <- heat.colors(5)
p <- ggplot(gap2, aes(x=gdpPercap, y=lifeExp)) +
  geom_point() +
  xlab("GDP per capita (US$, inflation-adjusted)") +
  ylab("Life expectancy (years)") +
  geom_point(data=gapminder[gapminder$country =="Kuwait",], col="gray") + 
  annotate("text", x=80000, y=70, label="Kuwait", col="gray")

## -------------------------------------------------------------------
wide.lm <- lm(lifeExp ~ poly(gdpPercap,5), data=gap2)
yilab <- c(0, 50, 100, 65, 85)
for (i in 2:5) {
    poly.lm <- lm(lifeExp ~ poly(gdpPercap,i), data=gap2)
    ft <- data.frame(x=gap2$gdpPercap, y=fitted(poly.lm))
    ft <- ft[order(ft$x),]
    p <- p + 
      geom_line(data=ft, aes(x=x,y=y), col=pal[i], lwd=2, alpha=0.8) + 
      annotate("text", x=60000, y=yilab[i], col="black", 
               label=sprintf("Polynomial degree %s", i))
}
gdp_focus <- c(10000, 25000, 40000)
p <- p + 
  geom_vline(xintercept=gdp_focus, col="gray") + 
  scale_x_continuous(breaks=c(0, gdp_focus, 60000, 80000, 100000))
p

## -------------------------------------------------------------------
inds <- rbind("quadratic"= c(1,1,1,0,0,0),
              "cubic"     =c(1,1,1,1,0,0),
              "quartic"   =c(1,1,1,1,1,0),
              "degree 5"  =c(1,1,1,1,1,1))
X <- newdata_to_X(list(gdpPercap=gdp_focus), wide.lm, intercept=TRUE)
rownames(X) <- gdp_focus
(ficres <- fic(wide.lm, inds=inds, focus=mean_normal, X=X))
summary(ficres)

## -------------------------------------------------------------------
median_normal<- function(par,X){
    qnorm(0.5, mean = as.numeric(X %*% par))
}
ficres <- fic(wide.lm, inds=inds, focus=median_normal, X=X)

## -------------------------------------------------------------------
q10_normal <- function(par, X, sigma){
    qnorm(0.1, mean = as.numeric(X %*% par), sd=sigma)
}
ficres <- fic(wide.lm, inds=inds, focus=q10_normal, X=X)

## -------------------------------------------------------------------
quantile_normal <- function(par, X, sigma, focus_p=0.5){
    qnorm(focus_p, mean = as.numeric(X %*% par), sd=sigma)
}

## -------------------------------------------------------------------
ficres <- fic(wide.lm, inds=inds, focus=quantile_normal,
               X=X[1,], focus_p=c(0.1,0.5,0.9))

## -------------------------------------------------------------------
focus_loglik <- function(par,X,sigma,Y){
    mu <- as.numeric(X %*% par)
    dnorm(Y,mu,sigma,log=TRUE)
}

## -------------------------------------------------------------------
wide.lm <- lm(mpg ~ am + wt + qsec + disp + hp, data=mtcars)
ncovs_wide <- length(coef(wide.lm)) - 1
inds0 <- c(1, rep(0, ncovs_wide))
inds <- all_inds(wide.lm, inds0)
X <- model.matrix(wide.lm)
Y <- model.response(model.frame(wide.lm))
ficres <- fic(wide.lm, inds=inds, focus=focus_loglik, X=X, Y=Y)

## -------------------------------------------------------------------
ficres <- ficres[ficres$vals=="ave",]
aics <- sapply(attr(ficres,"sub"), AIC)
qplot(ficres$rmse, aics,
      xlab="Root mean square error of log density estimate",
      ylab="AIC") 

