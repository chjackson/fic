* Check results vector includes everything useful, check formulae 

* return AIC, BIC if available for built in classes

* Testing against well understood examples, e.g. in book


* Model averaging weights


* More model classes

	- survival 
	- ??? 


* Analytic derivatives for particular focuses

  - interface for this implemented, needs a few more examples and testing
  
  - expected value in normal linear regression

* Vectorised FIC calculation.  I suggest to vectorise over different focuses rather than different submodels (e.g. different covariate categories)

  - This should allow averaged/weighted FIC over distribution of covariates. 

* Interface for comparing multiple models, e.g. all combinations

  - with some graphics 
  - inds argument as a matrix


* Cox models (needs additional theory?)


* Error handling for all functions



* Vignette + paper

	- Define clearly what model selection problems are supported. e.g. models that are nested in a common wide model? 
