##' Built-in focus functions and their derivatives
##' @name focus_fns
##'
##' @aliases prob_logistic prob_logistic_deriv
##'
##' @param par Vector of parameter estimates, including the intercept.
##'
##' @param X Vector or matrix of covariate values, including the intercept. This can either be a vector of length \eqn{p}, or a \eqn{n x p} matrix, where \eqn{p} is the number of covariate effects, and \eqn{n} is the number of alternative sets of covariate values at which the focus function is to be evaluated.
##'
##' @return \code{prob_logistic} returns the probability of the outcome in a logistic regression model, and \code{mean_normal} returns the mean outcome in a normal linear regression.   The \code{_deriv} functions return the vector of partial derivatives of the focus with respect to each parameter (or matrix, if there are multiple foci).
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
    t(X * p0 * (1 - p0))
}

##' @export
##' @rdname focus_fns
mean_normal <- function(par, X){
    X %*% par
}

##' @export
##' @rdname focus_fns
mean_normal_deriv <- function(par, X){
    t(X)
}

focus_fns <- 
    list(
        "prob_logistic" =
         list( 
             focus = prob_logistic,
             deriv = prob_logistic_deriv
         ),
       "mean_normal" =
         list( 
             focus = mean_normal,
             deriv = mean_normal_deriv
         )
       )

## Construct the arguments to supply to fic() from the built-in focus functions

get_focus <- function(focus, focus_deriv=NULL, par=NULL, X=NULL, ...){
    if (is.character(focus)) {
        if (!(focus %in% names(focus_fns)))
            stop(sprintf("focus function `%s` not found", focus))
        fi <- match(focus, names(focus_fns))
        args_extra <- list(...)
        focus <- focus_fns[[fi]]$focus
        deriv_args <- c(par=list(par), X=list(X), args_extra)
        focus_deriv <- do.call(focus_fns[[fi]]$deriv, deriv_args)
    }
    ## add unused X argument to function if it doesn't already have one
    if (!("X" %in% names(formals(focus))))
        formals(focus) <- c(formals(focus), alist(X=))
    list(focus=focus, focus_deriv=focus_deriv)
}
