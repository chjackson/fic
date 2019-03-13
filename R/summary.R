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
##' @export
summary.fic <- function(object, tidy=TRUE, adj=FALSE, ...) {
    if (adj)
        object$rmse <- object$rmse.adj
    minfic <- function(x){
        ind <- which.min(x$rmse)
        focus <- x$focus[ind]
        parnames <- parnames[inds[ind,]==1]
        if (tidy)
            data.frame(index=ind, focus=focus, pars=paste(parnames, collapse=","))
        else 
            list(index=ind, focus=focus, pars=parnames) 
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

## TODO test if focus not included 

print.summary.fic <- function(x){
    cat("Model with lowest RMSE by focus\n")
    print(x$min)
    cat("\nRange of focus estimates and RMSE over models\n")
    print(x$ranges)
}
