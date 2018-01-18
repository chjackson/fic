##' Built-in focus functions and their derivatives
##' @name focus_fns
##'
##' @aliases prob_logistic prob_logistic_deriv
##'
##' @param par Vector of parameter estimates
##'
##' @param X Vector or matrix of covariate values.  This can either be a vector of length \eqn{p}, or a \eqn{n x p} matrix, where \eqn{p} is the number of covariate effects, and \eqn{n} is the number of alternative sets of covariate values at which the focus function is to be evaluated.
NULL

##' @export
##' @rdname focus_fns
prob_logistic <- function(par, X){
    plogis(q = X %*% par)
}

##' @export
##' @rdname focus_fns
prob_logistic_deriv <- function(par, X){
    p0 <- as.vector(plogis(q = X %*% par))
    X * p0 * (1 - p0)
}

focus_fns <- 
  list("prob_logistic" =
         list( 
             focus = prob_logistic,
             deriv = prob_logistic_deriv
         )
       )

## Construct the arguments to supply to fic() from the built-in focus functions

get_focus <- function(focus, focus_deriv, par=NULL, ...){
    if (is.character(focus)) {
        if (!focus %in% names(focus_fns))
            stop(sprintf("focus function \"%s\" not found", focus))
        fi <- match(focus, names(focus_fns))
        args_extra <- list(...)
        focus <- focus_fns[[fi]]$focus
        deriv_args <- c(par=list(par), args_extra)
        focus_deriv <- do.call(focus_fns[[fi]]$deriv, deriv_args)
    } 
    list(focus=focus, focus_deriv=focus_deriv)
}
