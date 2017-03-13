myPlotFACETS <- function (x, emfit = NULL, clustered = FALSE, plot.type = c("none", "logR"), sname = NULL)
{
    def.par <- par(no.readonly = TRUE)
    plot.type <- match.arg(plot.type)
    if (plot.type == "none")
        layout(matrix(1:2, ncol = 1))
    if (plot.type == "em")
        layout(matrix(rep(1:4, c(9, 9, 6, 1)), ncol = 1))
    if (plot.type == "naive")
        layout(matrix(rep(1:4, c(9, 9, 6, 1)), ncol = 1))
    if (plot.type == "both")
        layout(matrix(rep(1:6, c(9, 9, 6, 1, 6, 1)), ncol = 1))
    par(mar = c(0.25, 3, 0.25, 1), mgp = c(1.75, 0.6, 0), oma = c(3,
        0, 1.25, 0))
    jseg <- x$jointseg
    if (missing(emfit)) {
        out <- x$out
        if (plot.type == "em" | plot.type == "both") {
            warning("emfit is missing; plot.type set to naive")
            plot.type <- "naive"
        }
    }
    else {
        out <- emfit$cncf
    }
    if (clustered) {
        cnlr.median <- out$cnlr.median.clust
        mafR <- out$mafR.clust
        mafR[is.na(mafR)] <- out$mafR[is.na(mafR)]
    }
    else {
        cnlr.median <- out$cnlr.median
        mafR <- out$mafR
    }
    mafR <- abs(mafR)
    chrcol <- 1 + rep(out$chr - 2 * floor(out$chr/2), out$num.mark)
    nn <- cumsum(table(jseg$chrom[is.finite(jseg$cnlr)]))
    segbdry <- cumsum(c(0, out$num.mark))
    segstart <- segbdry[-length(segbdry)]
    segend <- segbdry[-1]
    plot(jseg$cnlr[is.finite(jseg$cnlr)], pch = ".", cex = 0.5,
        col = c("grey", "lightblue", "azure4", "slateblue")[chrcol],
        ylab = "log-ratio", xaxt = "n", cex.axis=0.8, las=2)
    abline(h=0, col="lightgrey", lty=2)
    #abline(h = median(jseg$cnlr, na.rm = TRUE), col = "green2")
    #abline(h = x$dipLogR, col = "magenta4")
    segments(segstart, cnlr.median, segend, cnlr.median, lwd = 1.75,
        col = 2)
    if(plot.type!="logR"){
    plot(jseg$valor[is.finite(jseg$cnlr)], pch = ".", cex = 2.5,
        col = c("grey", "lightblue", "azure4", "slateblue")[chrcol],
        ylab = "log-odds-ratio", ylim = c(-4, 4), xaxt = "n", cex.axis=0.8, las=2)
    segments(segstart, sqrt(mafR), segend, sqrt(mafR), lwd = 1.75,
        col = 2)
    segments(segstart, -sqrt(mafR), segend, -sqrt(mafR), lwd = 1.75,
        col = 2)
    }
    chromlevels <- x$chromlevels
    if (is.null(chromlevels))
        chromlevels <- 1:length(nn)
    axis(labels = chromlevels, side = 1, at = (nn + c(0, nn[-length(nn)]))/2,
        cex.axis = 0.6)
    mtext(side = 1, line = 1.75, "Chromosome", cex = 1)
    if (!missing(sname))
        mtext(sname, side = 3, line = 0, outer = TRUE, cex = 0.8)
    par(def.par)
}
