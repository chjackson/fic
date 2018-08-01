##' Plot focused model comparison statistics
##'
##' Plot focused model comparison statistics.  What is plotted depends
##' on whether \code{\link{fic}} was called with a `sub` argument so that the
##' focus estimates are available in the object returned by \code{\link{fic}}.
##'
##' If the focus estimates are available, then the focus estimates are
##' plotted against the root MSE.  One plot is made for each covariate
##' value defining different focuses.
##'
##' If the focus estimates are unavailable, then the standard errors
##' of the focus estimate are plotted against the corresponding bias.
##' The plot points are shaded with darkness proportional to the RMSE,
##' with the point of maximum RMSE in black.
##'
##' @param x Output from \code{\link{fic}}.
##'
##' @param ylab y=axis label.
##'
##' @param xlab x-axis label.
##'
##' @param pch Plot point character, by default 19 (solid circle).
##'
##' @param mfrow Vector of two numbers giving the number of rows and number of columns respectively in the plot grid, if there are multiple focuses.
##' 
##' @param ... Other options to pass to \code{\link{plot}}.
##'
##' @export
plot.fic <- function(x, ylab=NULL, xlab=NULL, pch=19, mfrow=NULL, ...){
    if (is.array(x))
        x <- tidy.array(x, dim2=2, ord=c("vals","mods"))
    if (is.null(pch)) pch <- 19
    nfocus <- length(unique(x$vals))
    if (is.null(mfrow)) {
        mfrow <- grDevices::n2mfrow(nfocus)
    }
    par(mfrow=mfrow)

    for (i in seq(length=nfocus)){
        val <- unique(x$vals)[i]
        xi <- x[x$vals==val,,drop=FALSE]
        if (is.null(xi$focus)){
            prmse <- xi$rmse.adj / max(xi$rmse.adj)
            if (is.null(ylab)) ylab <- "SE"
            if (is.null(xlab)) xlab <- "Bias"
            plot(xi$bias.adj, xi$se, type="n", xlab=xlab, ylab=ylab, ...)
            abline(v=0.0, col="gray")
            points(xi$bias.adj, xi$se, pch=pch, col=gray(prmse))
        } else { 
            if (is.null(ylab)) ylab <- "Focus"
            if (is.null(xlab)) xlab <- "RMSE"
            plot(xi$rmse.adj, xi$focus, xlab=xlab, ylab=ylab, pch=pch, ...)
        }
        title(val)
    }
}

## TODO indicate wide and/or narrow models if stored in object
