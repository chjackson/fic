* Check results vector includes everything useful, check formulae 

* Testing against well understood examples, e.g. in book

* More model classes

   - survival
     - focus as extrapolated prob of survival.  document/suggest from flexsurv. short note in mdm?

	- illustrate for problems other than covariate selection, e.g. flexible survival models with shape parameter = 1
	- ???
    - illustrate for LMs: expected outcome at given cov val, particular cov
    eff, quantile at given cov val
	- other ones from book:
    	- polynomial regression in logistic
	    - skewed regression
		- poisson rate models
		- polynomial reg with het errors 
	- algebra for when the narrow model will be preferred 

* Cox models (Gerda has code for this)

* Analytic derivatives for particular focuses

  - interface for this implemented, needs a few more examples and testing
  - expected value in normal linear regression
  - though I doubt how useful this is

* Interface for comparing multiple models, e.g. all combinations

  - with some graphics 
  - inds argument as a matrix? 
  - or is clear example code enough?

* Illustrate averaged/weighted FIC

* Post-selection CIs.  Model averaging weights

* alternative losses.  prediction of an event

* Error handling for all functions

* testthat

* Vignette + paper


