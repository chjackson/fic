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
##'
##' @examples
##'
##' ## Model and focus from the main vignette 
##' wide.glm <- glm(low ~ lwtkg + age + smoke + ht + ui +
##'                 smokeage + smokeui, data=birthwt, family=binomial)
##' vals.smoke <-    c(1, 58.24, 22.95, 1, 0, 0, 22.95, 0)
##' vals.nonsmoke <- c(1, 59.50, 23.43, 0, 0, 0, 0, 0)
##' X <- rbind("Smokers" = vals.smoke, "Non-smokers" = vals.nonsmoke)
##' prob_logistic(coef(wide.glm), X=X)
##' prob_logistic_deriv(coef(wide.glm), X=X)
##'
##' @seealso \code{\link{fic}}
##' 
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
##'
##' @examples
##'
##' ## Mean mpg for a particular covariate category in the Motor Trend data
##' ## See the "fic" linear models vignette for more detail 
##' wide.lm <- lm(mpg ~ am + wt + qsec + disp + hp, data=mtcars)
##' cmeans <- colMeans(model.frame(wide.lm)[,c("wt","qsec","disp","hp")])
##' X <- rbind(
##'   "auto"   = c(intercept=1, am=0, cmeans),
##'   "manual" = c(intercept=1, am=1, cmeans)
##' )
##' mean_normal(coef(wide.lm), X)
##' mean_normal_deriv(coef(wide.lm), X)
##' 
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

get_focus <- function(focus, focus_deriv=NULL, par=NULL, ...){
    if (is.character(focus)) {
        if (!(focus %in% names(focus_fns)))
            stop(sprintf("focus function `%s` not found", focus))
        fi <- match(focus, names(focus_fns))
        args_extra <- list(...)
        focus <- focus_fns[[fi]]$focus
        deriv_args <- c(par=list(par), args_extra)
        focus_deriv <- do.call(focus_fns[[fi]]$deriv, deriv_args)
    }
    list(focus=focus, focus_deriv=focus_deriv)
}


check_focus <- function(focus, par, auxpar, ...)
{
    focus_args <- list(par=par)
    eargs <- list(...)
    nargs <- length(eargs)
    anames <- names(formals(focus))
    if (!(anames[1] == "par")) stop(sprintf("First argument of focus function named `%s`, this should be named `par`",anames[1]))
    focus_has_auxargs <- any(!(anames %in% c("par",names(eargs))))
    if (focus_has_auxargs)
        focus_args <- c(focus_args, auxpar)
    if (nargs > 0) { 
        for (i in seq_along(eargs)){
            eargs[[i]] <- check_focusarg(eargs[[i]], names(eargs)[i], length(par))
            ## on output these should all be matrices, nrow = nfocus, ncol = (typically) ncovs. no restriction on ncols
        }
        arglens <- sapply(eargs, nrow)
        nfocus <- max(arglens)
        ## Rownames of biggest argument used to label the different foci in the fic output
        wm <- which.max(arglens)
        if (is.null(rownames(eargs[[wm]])) | 
            any(duplicated(rownames(eargs[[wm]]))))
            rownames(eargs[[wm]]) <- eargs[[wm]][,1]
        for (i in seq_along(eargs)){
            if (arglens[i] == 1)
                eargs[[i]] <- eargs[[i]][rep(1,nfocus),,drop=FALSE]
            if (!(arglens[i] %in% c(1,nfocus)))
            {
                stop(sprintf("Number of focuses taken to be %s, but focus argument `%s` has %s elements/rows. This argument should either be a scalar or a vector/matrix with %s elements/rows", nfocus, names(eargs)[i], arglens[i], nfocus))
            }
        }
        fnames <- rownames(eargs[[wm]])
    } else {
        nfocus <- 1
        fnames <- NULL
        eargs <- NULL
    }
    focus_args <- c(focus_args, eargs)[anames]
    fn_try <- try(do.call(focus, focus_args))
    if (inherits(fn_try, "try-error")){
        argnames <- paste(paste0("`", names(formals(focus)[-1]), "`"), collapse=",")    
        stop(sprintf("Evaluating the focus function returned an error.  Check that the arguments %s to the focus function have been supplied and are valid, e.g. have the right dimensions", argnames))
    } 
    list(eargs=eargs, nfocus=nfocus, fnames=fnames)
}

check_focusarg <- function(arg, name, npar){
    if (!is.null(arg)){
        if (is.list(arg))
            arg <- as.matrix(as.data.frame(arg))
        if (!is.numeric(arg))
            stop(sprintf("Focus argument `%s` must be numeric", name))
        if (!(is.vector(arg) || is.matrix(arg))) 
            stop(sprintf("Focus argument `%s` must be a vector or a matrix", name))
        if (is.vector(arg)){
            if (name == "X")
                arg <- matrix(arg, nrow=1)
            else 
                arg <- matrix(arg, ncol=1)
        }
    }
    arg
}

