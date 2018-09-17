Hi Gerda,

Good news on the Cox regression.   After fixing a couple of bugs in my code, it now creates a FIC plot for the melanoma data (vignette p17) which looks just the same as figure 6.5 in your book. 

So I think my code is now doing what the code you used for the book is doing.

The only thing I don't understand is that the FICs produced by your code in this example are in a range between 0 and 1, rather than 0 to 30 in the book.   Is the code using a different definition of the FIC?  The focus values are identical I think, it's just the FIC values which are different.

If this can be explained simply, I'm happy to declare the Cox regression stuff finished and move on to another topic!

Chris





* rows of combs should be labelled by terms, not parameters

* error for wrong X dim in fic.coxph 

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
