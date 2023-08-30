#!/usr/bin/env Rscript

# detect script being run by snakemake
# if so, make a mock commandArgs function
if ('snakemake' %in% ls()) {
    logfile <- snakemake@log[[1]]
    con <- file(logfile, 'w')
    sink(con, type='output')
    sink(con, type='message')

    commandArgs <- function(...) {
        ret <- unlist(c(
            snakemake@input[1],         # table from make_roadmap_enrichment_table.R
            snakemake@params['ignore'],
            snakemake@output[1:3]      # output svg/pdf/csv ||no longer: /jpeg
        ))
        ret
    }
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 5) {
    cat("enrichment_table.csv is expected to have metadata columns in the order:\n")
    cat("  1. QUANTILES, 2. BINSIZE, 3-(n), metadata, (n+1)-end, esummary() rows.\n")
    cat("Any number of metadata columns is allowed; metadata column names cannot\n")
    cat("collide with esummary() column names.\n")
    cat("group_factors name metadata columns that are used to create plot columns.\n")
    cat("Every unique combination of group_factors will create a new column of plots.\n")
    cat("Other metadata will not affect columns.\n")
    cat("N.B. in these enrichment analyses, it is often better to group by 'datasource'\n")
    cat("rather than no grouping; this causes a title to be printed.\n")
    cat("ignore can either be the special string 'enrichment_grid_R_ignore_none' or\n")
    cat("a comma-separated list of key-value pairs of metadata to exclude from plotting.\n")
    cat("N.B., if BINSIZE is specified in ignore, then a blank box will be plotted to\n")
    cat("preserve structure with plots from other datasets.\n")
    stop("usage: make_roadmap_enrichment_grid.R enrichment_table.csv ignore out.svg out.pdf out.csv")
}

in.csv <- args[1]
ignore.string <- args[2]
out.svg <- args[3]
out.pdf <- args[4]
out.csv <- args[5]

if (file.exists(out.svg))
    stop(paste('output file', out.svg, 'already exists, please delete it first'))
if (file.exists(out.pdf))
    stop(paste('output file', out.pdf, 'already exists, please delete it first'))
if (file.exists(out.csv))
    stop(paste('output file', out.csv, 'already exists, please delete it first'))

library(data.table)
library(mutenrich)

max.plot.columns <- 100  # R's layout function maxes out at 200 columns

d <- fread(in.csv)
print(d)

# Get binsizes first so that if the user ignores them they still get plotted.
options(scipen=1000)
binsizes <- as.character(sort(as.integer(unique(d$BINSIZE))))
nbins <- length(binsizes)
print(binsizes)

# Get rid of the 'excluded' and 'outside' classes, which have non-integer quantile values
d <- d[!is.na(as.integer(d$quantile))]

if (ignore.string != 'enrichment_grid_R_ignore_none') {
    print(strsplit(ignore.string, ',')[[1]])
    for (i in strsplit(ignore.string, ',')[[1]]) {
        print(i)
        j <- strsplit(i, '=')[[1]]
        k <- j[1]
        v <- j[2]
        before <- nrow(d)
        d <- d[d[[k]] != v]
        after <- nrow(d)
        cat(paste0('Ignoring: ', before-after, 'rows with column=', k, ', value=', v, '\n'))
    }
} else {
    cat("Not ignoring any metadata columns.\n")
}

meta.cols <- colnames(d)[!(colnames(d) %in% c('QUANTILES', 'BINSIZE'))]
meta.cols <- meta.cols[(1:length(meta.cols)) < which(meta.cols=='quantile')]
# if "sigclass" is present (from sigenrich analysis), use it as the
# primary data grouping factor.
if ('sigclass' %in% meta.cols)
    meta.cols <- c('sigclass', meta.cols[-which(meta.cols=='sigclass')])
group.factors <- paste(meta.cols, collapse=',')
groups <- list(no.name=d)
if (length(group.factors) > 0) {
    group.factors <- strsplit(group.factors, ',')[[1]]
    # Determine the columns by group.factors
    cat('group.factors', group.factors, '\n')
    for (gf in group.factors) {
        if (!(gf %in% colnames(d)))
            stop(paste('group.factor', gf, 'is not a valid column name.\nAvailable column names:', colnames(d)))
    }
    groups <- split(d, d[, apply(.SD, 1, paste, collapse='\n'), .SDcols=group.factors])
}

cat("Got", length(groups), "groups (each = 1 plot column):\n")
print(names(groups))
ngroups <- min(max.plot.columns,length(groups))

# just make a blank panel with text 'txt' in the middle
textpane <- function(txt) {
    par(mar=c(0,0,0,0))
    plot(1, pch=NA, bty='n', xaxt='n', yaxt='n')
    legend('center', legend=txt, bty='n', cex=1.4)
}

suppressMessages(library(svglite))


# get height (in lines) of longest class name (will be in x-axis label)
em <- strheight('M', unit='inches')  # height of a capital M in inches
# 4 - y axis, 1 - right margin, 1.5 inches per panel (plus text label panel)
figwidth <- (4 + 1)*em + (1+ngroups)*1.5  # in inches

header.row <- max(strheight(names(groups), unit='inches'))*1.1 # 10% extra padding
# 2em margins top and bottom, 1.5 inches per panel + header row (column labels)
figheight <- nbins*1.5 + header.row + 4*em


# Loop over devices to save both pdf and svgs
devs <- list(pdf) #svglite, pdf)
outs <- c(out.pdf) #out.svg, out.pdf)
for (i in 1) { #1:2) {
    cat('making', outs[i], '\n')
    devs[[i]](file=outs[i], width=figwidth, height=figheight)

    layout(matrix(1:((nbins+1)*(ngroups+1)), nrow=(nbins+1)),
        heights=c(header.row,rep(1.5,nbins)))
    par(oma=c(2,2,2,2))
    for (bs in c('', binsizes)) {
        textpane(bs)
    }
    #for (group.name in names(groups)) {
    # by index now to facilitate splitting over multiple pages
    for (idx in 1:length(groups)) {
        group.name <- names(groups)[idx]

        # start new page
        if (idx %% max.plot.columns == 0) {
            layout(matrix(1:((nbins+1)*(ngroups+1)), nrow=(nbins+1)),
                heights=c(2,rep(4,nbins)))
            par(oma=c(2,2,2,2))
            for (bs in c('', binsizes)) {
                textpane(bs)
            }
        }
        textpane(strsplit(group.name, '\n')[[1]])
        for (bs in binsizes) {
            par(mar=c(1/2,3,0,0))
            g <- groups[[group.name]]
            g <- g[g$BINSIZE == bs]
            meta.cols <- colnames(g)[!(colnames(g) %in% c('QUANTILES', 'BINSIZE'))]
            meta.cols <- meta.cols[(1:length(meta.cols)) < which(meta.cols=='quantile')]
            signals <- split(g, g[, apply(.SD, 1, paste, collapse=' '), .SDcols=meta.cols])
            if (length(signals) == 0) {
                # Plot empty box placeholder
                plot(1, pch=NA, xaxt='n', yaxt='n')
            } else {
                for (i in seq_along(signals)) {
                    s <- signals[[i]]
                    s <- s[order(as.integer(s$quantile)),]
                    es <- data.frame(s)
                    rownames(es) <- s$quantile
                    bootstrap.ci <- if ('n.bootstraps' %in% colnames(es)) 0.95 else FALSE
                    plot.enrich(es=es,
                        bootstrap.ci=bootstrap.ci,
# svglite doesn't support multipage output
                        type='b',
                        show.asterisks=TRUE, cex.asterisk=2,
                        xaxt=ifelse(bs==tail(binsizes,1), 's', 'n'))
                        # not using uniform y axis anymore, so always have to plot it
                        #yaxt=ifelse(group.name==names(groups)[1], 's', 'n'))
                }

                abline(h=1, lty='dashed', col='grey')
            }
        }
    }
}

fwrite(d, file=out.csv)

if ('snakemake' %in% ls()) {
    sink()
}
