#!/usr/bin/env Rscript

# detect script being run by snakemake
# if so, make a mock commandArgs function
if ('snakemake' %in% ls()) {
    logfile <- snakemake@log[[1]]
    con <- file(logfile, 'w')
    sink(con, type='output')
    sink(con, type='message')

    commandArgs <- function(...) unlist(c(
        snakemake@input['oligo_snv_pta'],
        snakemake@input['oligo_snv_mda'],
        snakemake@output['svg'],
        snakemake@output['pdf'],
        snakemake@output['csv'],
        snakemake@params['signal_to_plot'],
        snakemake@params['features']
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 7) {
    cat('signal_to_plot: specified as a comma-separated key=value string. refers to the metadata in the asociated MANIFEST file.')
    cat('features: comma-separated list of features to plot along the x-axis. features are plotted in the order supplied.\n')
    stop("usage: suppfig2_panel_c.R oligo_pta_snv.enrich.summary.rda oligo_mda_snv.enrich.summary.rda out.svg out.pdf out.csv signal_to_plot features")
}

gpfile <- args[1]
gmfile <- args[2]
outsvg <- args[3]
outpdf <- args[4]
outcsv <- args[5]
signal_to_plot_string <- args[6]
features_string <- args[7]


if (file.exists(outsvg))
    stop(paste('output file', outsvg, 'already exists, please delete it first'))
if (file.exists(outpdf))
    stop(paste('output file', outpdf, 'already exists, please delete it first'))
if (file.exists(outcsv))
    stop(paste('output file', outcsv, 'already exists, please delete it first'))

suppressMessages(library(GenomicRanges))
suppressMessages(library(data.table))
suppressMessages(library(svglite))
suppressMessages(library(mutenrich))


load(gpfile)
gp <- es
load(gmfile)
gm <- es
# the last loaded emeta (which would be from this file, if it exists, is used)
# and it is assumed all emetas match.


if (!('list' %in% class(es))) {
    cat('ignoring signal_to_plot: enrichment files contain a single signal\n')
    signal_to_plot_string <- 'no signal selection'
    eslist <- list(
        `Oligo PTA SNVs`=gp,
        `Oligo MDA SNVs`=gm)

} else {
    signal.selections <- do.call(rbind, strsplit(strsplit(signal_to_plot_string, ',')[[1]], '='))
    colnames(signal.selections) <- c('key', 'value')
    cat('selecting signal via\n')
    print(signal.selections)
    emeta.idxs <- rep(TRUE, nrow(emeta))
    for (i in 1:nrow(signal.selections)) {
        cat('selecting:', signal.selections[i,1], '=', signal.selections[i,2],'\n')
        cat('passing before:', sum(emeta.idxs))
        emeta.idxs <- emeta.idxs & (emeta[,signal.selections[i,1]] == signal.selections[i,2])
        cat(' passing after:', sum(emeta.idxs), '\n')
    }
    if (sum(emeta.idxs) > 1) {
        warning(paste('only one signal will be plotted but', sum(emeta.idxs), 'matched the signal_to_plot criteria. using first matching signal'))
    }
    signal.to.keep <- which(emeta.idxs)[1]

    cat('keeping signal:', rownames(emeta)[signal.to.keep],'\n')

    eslist <- list(
        `Oligo PTA SNVs`=gp[[signal.to.keep]],
        `Oligo MDA SNVs`=gm[[signal.to.keep]])
}


feature.selections <- strsplit(features_string, ',')[[1]]
cat('selecting features for plotting (in order):', feature.selections, '\n')

eslist <- lapply(eslist, function(es) es[feature.selections,])


# ISSUES WITH GENERAL APPLICABILITY:
#   - need to specify par(mar), usually for x labels of classes
#   - need to specify class collapsing
#   - need to specify plot sizes

# all of the enrichment tables should have the same features, use this
# one for measurements.
x <- eslist[[1]]

print(rownames(x))
# get height (in lines) of longest class name (will be in x-axis label)
em <- strheight('M', unit='inches')  # height of a capital M in inches
xheight <- max(strwidth(rownames(x), unit='inches'))/em  # size in lines, for par(mar)

# 4 - y axis, 1 - right margin
figwidth <- 2*((4 + 1)*em + length(rownames(x))*em*2*1.5)  # in inches
# give 3.5 inches for the plot area, xheight for x-axis and 2 lines for title
# xheight*em = size in INCHES instead of lines
figheight <- 1*(3.5 + xheight*em + 2*em)  # 1* - one plot high. maybe in the future allow multiple signals in one plot by stacking plots vertically


# turn enrichment as a fraction (e.g., 0.92 or 1.22 representing 8% and
# 22% depletions and enrichments, respectively) into a positive or negative
# percent
scale.enr <- function(x) 100 * (x - 1)

# point.ests - list of bar heights
# bars - list of error bar lower and upper bounds (columns 1 and 2 of a matrix each)
plot.one <- function(point.ests, bars, ylim, cols=1:length(point.ests), ...) {
    points <- scale.enr(do.call(rbind, point.ests))
    # bp has locations of each bar for plotting error bars
    bp <- barplot(points, names.arg=names(point.ests[[1]]), beside=TRUE,
        las=3, col=cols, border=NA, ylim=ylim, ylab='Enrichment or depletion (%)', ...)
    abline(h=0)

    for (i in 1:nrow(bp)) {
        xloc <- bp[i,]
        arrows(x0=xloc, y0=scale.enr(bars[[i]][,1]), x1 = xloc, y1=scale.enr(bars[[i]][,2]),
            angle = 90, code = 3, length = 0.02, lwd = 1, col=cols[[i]])
    }
}

# Loop over devices to save both pdf and svgs
devs <- list(svglite, pdf)
outs <- c(outsvg, outpdf)
for (i in 1:2) {
    devs[[i]](file=outs[i], width=figwidth, height=figheight)

    par(mar=c(xheight, 4, 3, 1))
    layout(t(1:2))

    ylim <- range(c(get.ylim(es=eslist[[1]], shift=1, scale=100, use.boot=T),
                    get.ylim(es=eslist[[2]], shift=1, scale=100, use.boot=T)))
    plot.one(point.ests=list(setNames(eslist[[1]]$enr, rownames(eslist[[1]])), eslist[[2]]$enr),
        bars=lapply(eslist[1:2], function(es) es[, c("enr.boot.0.95.lb", "enr.boot.0.95.ub")]),
        cols=c('red', 'purple'),
        #ylim=ylim, family='Arial', main=paste0('OL sSNVs\n', signal_to_plot_string))
        ylim=ylim, main=paste0('OL sSNVs\n', signal_to_plot_string))

    plot(1, pch=NA, xaxt='n', yaxt='n', xlab='', ylab='', main='', bty='n')
    legend("top", fill=c('red', 'purple'), legend=c("PTA OLs", "MDA OLs"), bty='n')
    dev.off()
}


# Write table of statistics
stat.table <- do.call(rbind, lapply(1:length(eslist), function(i) {
    es <- eslist[[i]]
    label <- names(eslist)[i]
    cbind(label, class=rownames(es), es)
}))
fwrite(stat.table, file=outcsv)

if ('snakemake' %in% ls()) {
    sink()
}
