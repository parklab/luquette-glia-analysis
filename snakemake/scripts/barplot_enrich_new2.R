#!/usr/bin/env Rscript

# detect script being run by snakemake
# if so, make a mock commandArgs function
if ('snakemake' %in% ls()) {
    logfile <- snakemake@log[[1]]
    con <- file(logfile, 'w')
    sink(con, type='output')
    sink(con, type='message')

    commandArgs <- function(...) unlist(c(
        snakemake@output['svg'],
        snakemake@output['pdf'],
        snakemake@params['colors'],
        snakemake@input['enrich_tables']
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 3) {
    cat("Assumes a consolidate_*.R script has been previously run.\n")
    cat("signal to plot should have already been selected and features both selected and ordered.\n")
    cat("All consolidated enrichment tables should have the same set of features in the same order.\n")
    cat("`colors` is a comma-separated string of R-recognized color strings. Order must match the consolidated_summary csv order\n")
    stop("usage: barplot_enrich_new2.R out.svg out.pdf colors consolidated_summary1.csv [ consolidated_summary2.csv ... consolidated_summaryN.csv ]")
}

outsvg <- args[1]
outpdf <- args[2]
colors <- unlist(strsplit(args[3], split=','))
print(colors)
in.csvs <- args[-(1:3)]

if (file.exists(outsvg))
    stop(paste('output file', outsvg, 'already exists, please delete it first'))
if (file.exists(outpdf))
    stop(paste('output file', outpdf, 'already exists, please delete it first'))

suppressMessages(library(data.table))
suppressMessages(library(svglite))
suppressMessages(library(mutenrich))

# turn enrichment as a fraction (e.g., 0.92 or 1.22 representing 8% and
# 22% depletions and enrichments, respectively) into a positive or negative
# percent
scale.enr <- function(x) 100 * (x - 1)

eslist <- lapply(in.csvs, function(infile) { 
    ret <- fread(infile)
    ret$enr <- scale.enr(ret$enr)
    ret$enr.boot.0.95.lb <- scale.enr(ret$enr.boot.0.95.lb)
    ret$enr.boot.0.95.ub <- scale.enr(ret$enr.boot.0.95.ub)
    ret
})
names(eslist) <- sapply(eslist, function(es) es$group.tag[1])

# all of the enrichment tables should have the same features, use this
# one for measurements.
xaxis.labels <- eslist[[1]]$class

print(xaxis.labels)
# get height (in lines) of longest class name (will be in x-axis label)
emh <- strheight('M', unit='inches')  # height of a capital M in inches
emw <- strwidth('M', unit='inches')  # height of a capital M in inches

# Try to estimate string width, when turned 90 degrees, in lines. R's
# strheight/strwidth do not seem to calcualte this case.
# factor of 0.8 - most characters are not as wide as capital M, try 85% as an average.
xheight <- max(strwidth(xaxis.labels, unit='inches'))/emw  * 0.85 # size in lines, for par(mar)
print(xheight)
xheight.in <- max(strwidth(xaxis.labels, unit='inches')) #/emw  * 0.8 # size in lines, for par(mar)
print(xheight.in)

# 4 - y axis, 1 - right margin
barplot.width <- (4 + 1)*emw + length(xaxis.labels)*emw*length(eslist)*1.5  # in inches
legend.width <- max(strwidth(names(eslist), unit='inches')) + 6*emw
figwidth <- barplot.width + legend.width

# give 3.5 inches for the plot area, xheight for x-axis and 2 lines for title
# 0*2*em - no longer plotting a title
# xheight*em = size in INCHES instead of lines
#figheight <- 1*(3.5 + xheight*emh + 0*2*emh)  # 1* - one plot high. maybe in the future allow multiple signals in one plot by stacking plots vertically
figheight <- 3.5 + xheight.in



# point.ests - list of bar heights
# bars - list of error bar lower and upper bounds (columns 1 and 2 of a matrix each)
plot.one <- function(point.ests, bars, ylim, cols=1:length(point.ests), ...) {
    points <- do.call(rbind, point.ests)
    # bp has locations of each bar for plotting error bars
    bp <- barplot(points, names.arg=names(point.ests[[1]]), beside=TRUE,
        las=3, col=cols, border=NA, ylim=ylim, ylab='Enrichment or depletion (%)', ...)
    abline(h=0)

    for (i in 1:nrow(bp)) {
        xloc <- bp[i,]
        arrows(x0=xloc, y0=bars[[i]][[1]], x1 = xloc, y1=bars[[i]][[2]],
            angle = 90, code = 3, length = 0.02, lwd = 1, col=cols[[i]])
    }
}

# Loop over devices to save both pdf and svgs
devs <- list(svglite, pdf)
outs <- c(outsvg, outpdf)
for (i in 1:2) {
    devs[[i]](file=outs[i], width=figwidth, height=figheight)

    layout(matrix(1:2, nrow=1), width=c(barplot.width, legend.width))

    par(mar=c(xheight, 4, 1, 1))

    ylim <- range(pretty(unlist(lapply(eslist, function(es) c(es$enr, es$enr.boot.0.95.lb, es$enr.boot.0.95.ub)))))
        print(lapply(eslist, function(es) es[, c("enr.boot.0.95.lb", "enr.boot.0.95.ub")]))
    plot.one(point.ests=lapply(eslist, function(es) setNames(es$enr, es$class)),
        bars=lapply(eslist, function(es) es[, c("enr.boot.0.95.lb", "enr.boot.0.95.ub")]),
        ylim=ylim, cols=colors)

    par(mar=c(1,1,1,1))
    plot(0, pch=NA, yaxt='n', xaxt='n', xlab='n', ylab='n', bty='n')
    legend('center', fill=colors, legend=names(eslist), bty='n')

    dev.off()
}

if ('snakemake' %in% ls()) {
    sink()
}
