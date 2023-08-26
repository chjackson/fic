fic
===

The development repository for the [fic](http://cran.r-project.org/package=fic) R package for the Focused Information Criterion and related methods for model comparison. 

The `fic` package compares how well different models estimate a quantity of interest (the "focus") so that different models may be preferred for different purposes.

Comparisons within any class of models fitted by maximum likelihood are supported, with shortcuts for commonly-used classes such as generalised linear models and parametric survival models.

The methods originate from [Claeskens and Hjort (2003)](https://amstat.tandfonline.com/doi/abs/10.1198/016214503000000819#.XLiKGy2ZOu4) and [Claeskens and Hjort (2008)](https://feb.kuleuven.be/public/u0043181/modelselection/)

## Installation (stable CRAN version)
```r
install.packages("fic")
```

## Installation (development version)

```r
install.packages("devtools") # if devtools not already installed
library(devtools)
install_github("chjackson/fic")
 ```

See the package vignettes (work in progress) for an introduction:

[Different models for different purposes: focused model comparison in R](https://chjackson.github.io/fic/inst/doc/fic.pdf)

[Linear models example](https://chjackson.github.io/fic/inst/doc/linear.pdf)

[Skew-normal example](https://chjackson.github.io/fic/inst/doc/skewnormal.pdf)

[Survival analysis example](https://chjackson.github.io/fic/inst/doc/survival.pdf)

[Multi-state models example](https://chjackson.github.io/fic/inst/doc/multistate.pdf)

[Bootstrap and alternative losses](https://chjackson.github.io/fic/inst/doc/loss.pdf)

Source code [GitHub repository](https://github.com/chjackson/fic)

<!-- badges: start -->
[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![R-CMD-check](https://github.com/chjackson/fic/workflows/R-CMD-check/badge.svg)](https://github.com/chjackson/fic/actions)
[![Codecov test coverage](https://codecov.io/gh/chjackson/fic/branch/master/graph/badge.svg)](https://app.codecov.io/gh/chjackson/fic?branch=master)
<!-- badges: end -->
