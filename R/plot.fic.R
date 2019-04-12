##' Plot focused model comparison statistics: base graphics method
##'
##' Plot focused model comparison statistics: base graphics method
##'
##' If the focus estimates are available, then the focus estimates are
##' plotted against the root MSE.  One plot is made for each covariate
##' value defining different focuses.  If the wide model estimate is
##' available, this is illustrated as a solid line on the plot, and if
##' the narrow model estimate is available, this is shown as a dashed
##' line.
##'
##' If the focus estimates are unavailable, then the standard errors
##' of the focus estimate are plotted against the corresponding bias.
##' The plot points are shaded with darkness proportional to the RMSE,
##' with the point of maximum RMSE in black.
##'
##' The \pkg{ggplot2}-based plot method, \code{\link{ggplot_fic}}, is
##' slightly nicer. 
##'
##' @param x Output from \code{\link{fic}}.
##'
##' @param ci Plot interval estimates? (\code{TRUE} or \code{FALSE}).  These are calculated as plus / minus twice the standard error of the submodel focus under the wide model.  These are rough estimates of uncertainty intended to illustrate the bias-variance tradeoff, and exclude any uncertainty associated with the choice between models.
##' 
##' @param xlab x-axis label.
##'
##' @param ylab y-axis label.
##'
##' @param xlim x-axis limits (pair of numbers)
##'
##' @param ylim y-axis limits
##'
##' @param pch Plot point character, by default 19 (solid circle).
##'
##' @param mfrow Vector of two numbers giving the number of rows and number of columns respectively in the plot grid, if there are multiple focuses.
##' 
##' @param ... Other options to pass to \code{\link{plot}}.
##'
##' @seealso \code{\link{ggplot_fic}}, \code{\link{summary.fic}}
##'
##' @import graphics grDevices
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
##' plot(ficres)
##'
##' @export
plot.fic <- function(x, ci=TRUE, xlab=NULL, ylab=NULL, xlim=NULL, ylim=NULL, pch=19, mfrow=NULL, ...){
    if (is.array(x))
        x <- tidy.array(x, dim2=2, ord=c("vals","mods"))
    if (is.null(pch)) pch <- 19
    nfocus <- length(unique(x$vals))
    if (is.null(mfrow)) {
        mfrow <- rev(grDevices::n2mfrow(nfocus))
    }
    par(mfrow=mfrow)

    if (!is.null(x$focus)){
        x$l95 = x$focus - qnorm(0.975)*x$se
        x$u95 = x$focus + qnorm(0.975)*x$se
    }
    iwide <- attr(x, "iwide")
    inarr <- attr(x, "inarr")
    col <- scales::hue_pal()(3)[3]

    for (i in seq(length=nfocus)){
        val <- unique(x$vals)[i]
        xi <- x[x$vals==val,,drop=FALSE]
        if (is.null(xi$focus)){
            prmse <- xi$rmse.adj / max(xi$rmse.adj)
            if (is.null(ylab)) ylab <- "SE"
            if (is.null(xlab)) xlab <- "Bias"
            plot(xi$bias.adj, xi$se, type="n",
                 xlim=xlim, ylim=ylim,
                 xlab=xlab, ylab=ylab, ...)
            abline(v=0.0, col="gray")
            points(xi$bias.adj, xi$se, pch=pch, col=gray(prmse))
        } else { 
            if (is.null(xlab)) xlab <- "Focus"
            if (is.null(ylab)) ylab <- "RMSE"
            if (is.null(xlim)) xlim <- range(c(x$l95, x$u95))
            if (is.null(ylim)) ylim <- range(x$rmse.adj)
            plot(xi$focus, xi$rmse.adj, type="n",
                 xlim=xlim, ylim=ylim,
                 xlab=xlab, ylab=ylab, ...)
            if (!is.null(iwide))
                abline(v=xi$focus[iwide], col=col)
            if (!is.null(inarr))
                abline(v=xi$focus[inarr], col=col, lty=2)
            points(xi$focus, xi$rmse.adj, pch=pch, ...)
            if (ci) 
                segments(xi$l95, xi$rmse.adj, xi$u95, xi$rmse.adj)
            text(min(x$l95), xi$rmse.adj, labels=xi$mods, col="gray60", cex=0.7, pos=4)
        }
        title(val)
    }
}


##' Plot focused model comparison statistics: ggplot2 method
##'
##' This only works if the focus estimates are available. The focus estimates are
##' plotted against the root MSE.  One plot is made for each covariate
##' value defining different focuses.  If the wide model estimate is
##' available, this is illustrated as a solid line on the plot, and if
##' the narrow model estimate is available, this is showm as a dashed
##' line.
##'
##' @inheritParams plot.fic
##'
##' @param legend Show a legend, identifying
##'
##' a) the line types for the wide and narrow models
##'
##' b) the names of the terms of the wide model.  This is used when
##' the \code{inds} object supplied to \code{\link{fic}} was constructed by
##' \code{\link{all_inds}}, so has row names made out of a string
##' of 0s and 1s that identify the terms included in the submodel.
##' These strings are plotted as text labels against the estimate for
##' each submodel.  The legend identifies which 0s and 1s correspond
##' to which model terms.
##'
##' @importFrom scales hue_pal
##' 
##' @import ggplot2
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
##' @seealso plot.fic, summary.fic
##' 
##' @export
ggplot_fic <- function(x, ci=TRUE, legend=TRUE, ylab=NULL, xlab=NULL, xlim=NULL, ylim=NULL){
    if (is.null(x$focus))
        stop("No focus estimates found. `fic` should be run with the `sub` argument")
    if (is.null(ylab)) ylab <- "RMSE"
    if (is.null(xlab)) xlab <- "Focus"
    x$l95 = x$focus - qnorm(0.975)*x$se
    x$u95 = x$focus + qnorm(0.975)*x$se
    iwide <- attr(x, "iwide")
    inarr <- attr(x, "inarr")
    col <- scales::hue_pal()(3)[3]
    ps <- ggplot(data=x, aes_string(x='focus', y='rmse.adj'))
    if (!is.null(iwide)){
        indwide <- tapply(1:nrow(x), x$vals, function(x)x[iwide])
        xwide <- x[indwide,,drop=FALSE]
        xwide$Model <- "Wide"
    }
    if (!is.null(inarr)){
        indnarr <- tapply(1:nrow(x), x$vals, function(x)x[inarr])
        xnarr <- x[indnarr,,drop=FALSE]
        xnarr$Model <- "Narrow"
    }
    xmod <- rbind(xnarr, xwide)
    xmod$Model <- factor(xmod$Model, levels=c("Wide","Narrow"))
    ps <- ps + geom_vline(data=xmod, aes_string(xintercept = 'focus', lty='Model'), col=col)
    ps <- ps + 
        facet_grid(.~vals) + 
        geom_point(aes(color="black")) + # to force legend
        scale_colour_identity(guide="legend") + # allow "black" in aes instead of variable
        xlab(xlab) +
        ylab(ylab)
    if (!is.null(xlim))
        ps <- ps + xlim(xlim)
    if (ci)
        ps <- ps +
            geom_segment(aes_string(x='l95', xend='u95', yend='rmse.adj'))
    ps <- ps +
      geom_text(aes_string(x=-Inf, label='mods', hjust=0), col="gray60", size=2.5)
    pn <- attr(x, "termnames")
    if (!is.null(pn) && legend){
        plab <- paste0("Terms\n", paste(paste(seq(along.with=pn), pn, sep=": "), collapse="\n"))
        ## Hack a legend to display parameter names 
        ps <- ps +
            guides(color = guide_legend(title=plab, # use par names as legend title
                                        label=FALSE, # remove key label 
                                        override.aes = list(alpha = 0))) + # remove key background
            theme(legend.key=element_blank()) # remove key point symbols
    } else ps <- ps + theme(legend.position="none")
    ps
}
