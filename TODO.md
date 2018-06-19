## Definitely

* More model classes and built-in focuses, to include at least 

	- survival vignette 

	- linear model vignette

	- skewed regression vignette

* Cox models

* Illustrate some algebra/intuitions for when the narrow model will be preferred? 

* narrow model pars in the middle - test, document 

* Handle regression intercepts more nicely?

* Illustrate averaged/weighted FIC
  - why can't we use weighted average of FIC - is there a bias which
    is corrected by doing the algebra on p180-181? 

* Alternative losses, e.g. error rate for prediction of an event. 
  - Works well empirically to compute replicate focuses using MVN sample from wide model, then use these to calculate MSE of submodel focus MLE. 
  - Assume this is equivalent to procedure on pp152-3 on book
  - Can substitute any other loss for MSE 
  - Todo discuss with Gerda

* Error handling for all functions
  - can narrow model parameters be in the middle? 
  - error if submodel doesn't match inds 

* testthat tests

* Vignettes and JSS paper


## Possibly
  
* High-dimensional regression 
 - "They present one method for when you can fit the submodel but not the wide model, and the criterion is computed using estimates from the submodel.  Another method is for when both submodel and wide model are high dimensional, and it uses penalisation."

* Post-selection CIs.  Model averaging weights
  - bootstrap solution?   bootstrap from the data, do FIC selection each time, then take empirical CI of focus estimate 
  - recommend in paper 
  
