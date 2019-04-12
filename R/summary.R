get_parnames <- function(par, inds){
    if (!is.null(names(par)))
        res <- names(par)
    else if (!is.null(colnames(inds)))
        res <- colnames(inds)
    else res <- NULL
    res
}    

##' Summarise focused model comparison results
##'
##' @param object Object returned by \code{\link{fic}} representing focused model comparison statistics for a range of models, and potentially also multiple focus quantities.
##'
##' @param tidy If \code{TRUE} (the default) then the results describing the optimal model (per focus) are returned as a data frame, with the names of the parameters in the optimal model collapsed into a single string.  If \code{FALSE}, the results are returned as a list, including a vector of parameter names.
##'
##' @param adj The optimal model is the one with the lowest root mean square error (RMSE). If \code{adj=TRUE} the RMSE is based on the adjusted bias estimator.  Otherwise the standard estimator is used. 
##'
##' @param ... Other arguments, currently unused. 
##'
##' @return A list of two components, one for the optimal model per focus, and one for the range of focus and RMSE estimates over models.
##'
##' @seealso \code{\link{ggplot_fic}}, \code{\link{plot.fic}} for a more detailed visual representation of the focused comparison
##'
##' @examples
##'
##' ## Example from the main vignette, see there for more details
##' 
##' wide.glm <- glm(low ~ lwtkg + age + smoke + ht + ui + smokeage + smokeui,
##'                 data=birthwt, family=binomial)
##' vals.smoke <-    c(1, 58.24, 22.95, 1, 0, 0, 22.95, 0)
##' vals.nonsmoke <- c(1, 59.50, 23.43, 0, 0, 0, 0, 0)
##' X <- rbind("Smokers" = vals.smoke, "Non-smokers" = vals.nonsmoke)
##' inds0 <- c(1,1,0,0,0,0,0,0)
##' combs <- all_inds(wide.glm, inds0)
##' ficres <- fic(wide = wide.glm, inds = combs, inds0 = inds0,
##'               focus = prob_logistic, X = X)
##' ggplot_fic(ficres)
##' summary(ficres)
##'
##' @export
summary.fic <- function(object, tidy=TRUE, adj=FALSE, ...) {
    if (adj)
        object$rmse <- object$rmse.adj
    if (is.null(object$focus)) {
        object$focus <- NA
        warning("Focus values not available")
    }
    minfic <- function(x){
        ind <- which.min(x$rmse)
        focus <- x$focus[ind]
        parnames <- parnames[inds[ind,]==1]
        if (tidy)
            res <- data.frame(index=ind, pars=paste(parnames, collapse=","))
        else 
            res <- list(index=ind, pars=parnames)
        res$focus <- focus 
        res
    }
    fsplit <- split(object, object$vals)
    parnames <- attr(object, "parnames")
    inds <- attr(object, "inds")
    resmin <- lapply(fsplit, minfic)
    if (tidy) resmin <- do.call("rbind", resmin)
    frange <- as.data.frame(t(sapply(fsplit, function(x) range(x$focus))))
    rrange <- as.data.frame(t(sapply(fsplit, function(x) range(x$rmse, na.rm=TRUE))))
    ranges <- cbind(frange, rrange)
    names(ranges) <- c("min(focus)", "max(focus)", "min(RMSE)", "max(RMSE)")
    res <- list(min=resmin, ranges=ranges)
    class(res) <- "summary.fic"
    res
}

print.summary.fic <- function(x){
    cat("Model with lowest RMSE by focus\n")
    print(x$min)
    cat("\nRange of focus estimates and RMSE over models\n")
    print(x$ranges)
}
