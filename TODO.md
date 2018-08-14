## Definitely

* More model classes and built-in focuses in vignette(s)

	- survival

	- linear models

	- skewed regression

	- multi-state models

* Cox models.  Code written but don't currently trust it

* Alternative losses, e.g. error rate for prediction of an event. 
  - Works well empirically to compute replicate focuses using MVN sample from wide model ests, then use these to calculate MSE of submodel focus MLE. 
  - Assume this is equivalent to procedure on pp152-3 on book
  - Can substitute any other loss for MSE 
  - Todo discuss with Gerda

* Illustrate some algebra/intuitions for when the narrow model will be preferred? 

* Make sure all help pages complete

* Error handling for all functions
  - use imagination here, go through all args

* testthat tests

* Vignettes and JSS paper


## Possibly
  
* High-dimensional regression 
 - "They present one method for when you can fit the submodel but not the wide model, and the criterion is computed using estimates from the submodel.  Another method is for when both submodel and wide model are high dimensional, and it uses penalisation."

* Post-selection CIs.  Model averaging weights
  - bootstrap solution?   bootstrap from the data, do FIC selection each time, then take empirical CI of focus estimate 
  - recommend in paper 
