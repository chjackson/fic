* Special case in fic for all parameters excluded

* Refine what is included in results vector

	- fic as root MSE
	- different kinds of bias/variance/fic estimators
	- estimated focus quantity under each submodel
	- AIC, BIC? would need to refit the submodel.
	- model averaging weights


* Testing against well understood examples, e.g. in book


* More model classes

	- survival 
	- ??? 


* Analytic derivatives for particular focuses

  - prob of outcome in logistic regression

  - expected value in normal linear regression

  - interface for this - character string, rather than name of a function?

  - should still be able to pass arguments, e.g. covariate values, to the focus
  

* Interface for comparing multiple models, e.g. all combinations

  - with some graphics 


* Averaged/weighted FIC over distribution of covariates


* Cox models (needs additional theory?)


* Error handling for all functions



* Vignette + paper

	- Define clearly what model selection problems are supported. e.g. models that are nested in a common wide model? 
