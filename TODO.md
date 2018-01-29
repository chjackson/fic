* Check results vector includes everything useful, check formulae 

* Testing against well understood examples, e.g. in book

* Utilities for forming "inds" matrices for all combinations of models, e.g. regression subsets. 

* More model classes

   - survival
     - focus as expected survival over horizon, as in HE mods. short note in mdm?

    - illustrate for LMs: expected outcome at given cov val, particular cov
    eff, quantile at given cov val

	- other ones from book:
    	- polynomial regression in logistic
	    - skewed regression
		- poisson rate models
		- polynomial reg with het errors 

	- ???

* Illustrate some algebra/intuitions for when the narrow model will be preferred? 

* Cox models (Gerda has code for this)

* Analytic derivatives for particular focuses

  - interface for this implemented, needs a few more examples and testing
  - expected value in normal linear regression
  - though I doubt how useful this is

* Illustrate averaged/weighted FIC

* Post-selection CIs.  Model averaging weights

* Alternative losses.  Prediction of an event

* Error handling for all functions

* testthat tests

* Vignette + paper
