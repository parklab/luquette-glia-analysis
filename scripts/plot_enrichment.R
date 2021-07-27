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
            snakemake@output[1:3],      # output svg/pdf/csv ||no longer: /jpeg
            snakemake@params['xlab']
        ))
        if ('xgroup' %in% names(snakemake@params))
            ret <- c(ret, unlist(snakemake@params['xgroup']))
        if ('ygroup' %in% names(snakemake@params)) {
            if (!('xgroup' %in% names(snakemake@params)))
                stop('if ygroup is specified, xgroup must also be specified')
            ret <- c(ret, unlist(snakemake@params['ygroup']))
        }
        if ('highlight' %in% names(snakemake@params)) {
            if (!('ygroup' %in% names(snakemake@params)))
                stop('if highlight is specified, ygroup must also be specified')
            ret <- c(ret, unlist(snakemake@params['highlight']))
        }
        ret
    }
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}

args <- commandArgs(trailingOnly=TRUE)
if (all(length(args) != 6:9)) {
    cat('This script is just a copy of plot_enrichment_grid.R with bins and quantiles removed and outside/excluded regions retained\n')
    stop("usage: make_roadmap_enrichment_grid.R enrichment_table.csv ignore out.svg out.pdf out.csv xlabfactor1,...,xlabfactorN [ x_group_factor1,...,x_group_factorN ] [ y_group_factor1,...,y_group_factorN ] [highlight]")
}

in.csv <- args[1]
ignore.string <- args[2]
out.svg <- args[3]
out.pdf <- args[4]
out.csv <- args[5]
xlab.factors <- args[6]
xgroup.factors <- c()
if (length(args) >= 7)
    xgroup.factors <- args[7]
ygroup.factors <- c()
if (length(args) >= 8)
    ygroup.factors <- args[8]
highlight <- c()
if (length(args) >= 9)
    highlight <- args[9]

if (file.exists(out.svg))
    stop(paste('output file', out.svg, 'already exists, please delete it first'))
if (file.exists(out.pdf))
    stop(paste('output file', out.pdf, 'already exists, please delete it first'))
if (file.exists(out.csv))
    stop(paste('output file', out.csv, 'already exists, please delete it first'))

library(data.table)

d <- fread(in.csv)

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

if (length(highlight) > 0) {
    highlight.key <- strsplit(highlight, '=')[[1]][1]
    highlight.value <- strsplit(highlight, '=')[[1]][2]
    if (!(highlight.key %in% colnames(d)))
        stop(paste('requested highlight metadata column', highlight.key, 'is not in the table. Did you ignore it?'))

    cat(paste('Highlighting', sum(d[[highlight.key]] == highlight.value),
        'instances of', highlight.key, '=', highlight.value, '\n'))
}


xlab.factors <- strsplit(xlab.factors, ',')[[1]]
for (gf in xlab.factors) {
    if (!(gf %in% colnames(d)))
        stop(paste('xlab.factor', gf, 'is not a valid column name.\nAvailable column names:', colnames(d)))
}
d[, x_label_factor_class := apply(.SD, 1, paste, collapse='\n'),
    .SDcols=xlab.factors]
xlabs <- unique(d[['x_label_factor_class']])
nxlabs <- length(xlabs)
cat("Got", nxlabs, "unique x labels (each = 1 point per plot):\n")
print(xlabs)


if (length(xgroup.factors) > 0) {
    xgroup.factors <- strsplit(xgroup.factors, ',')[[1]]
    for (gf in xgroup.factors) {
        if (!(gf %in% colnames(d)))
            stop(paste('xgroup.factor', gf, 'is not a valid column name.\nAvailable column names:', colnames(d)))
    }
    d[, x_group_factor_class := apply(.SD, 1, paste, collapse='\n'),
        .SDcols=xgroup.factors]
} else {
    d[, x_group_factor_class := '']
}
xgroups <- unique(d[['x_group_factor_class']])
nxgroups <- length(xgroups)
cat("Got", nxgroups, "x groups (each = 1 plot row):\n")
print(xgroups)


if (length(ygroup.factors) > 0) {
    for (gf in ygroup.factors) {
        if (!(gf %in% colnames(d)))
            stop(paste('ygroup.factor', gf, 'is not a valid column name.\nAvailable column names:', colnames(d)))
    }
    d[, y_group_factor_class := apply(.SD, 1, paste, collapse='\n'),
        .SDcols=ygroup.factors]
} else {
    d[, y_group_factor_class := '']
}
ygroups <- unique(d[['y_group_factor_class']])
nygroups <- length(ygroups)
cat("Got", nygroups, "y groups (each = 1 plot column):\n")
print(ygroups)

print(d)


# just make a blank panel with text 'txt' in the middle
textpane <- function(txt, legend.position='center') {
    par(mar=c(0,0,0,0))
    plot(1, pch=NA, bty='n', xaxt='n', yaxt='n')
    legend(legend.position, legend=txt, bty='n', cex=1.4)
}

suppressMessages(library(extrafont))
suppressMessages(library(svglite))
suppressMessages(library(mutenrich))


if (!("Arial" %in% fonts()))
    stop("Arial font not detected; did you load extrafonts and run font_import() with the appropriate path?")

# get height (in lines) of longest class name (will be in x-axis label)
em <- strheight('M', unit='inches')  # height of a capital M in inches
# 4 - y axis, 1 - right margin, 4 - outer margins
# 1.5em per point per y panel (plus text label panel)
figwidth <- (4 + 1 + 4)*em + max(strwidth(xgroups, unit='inches')) + nygroups*em*nxlabs*1.5  # in inches

# x labels will be rotated, so use width
emw <- strwidth('M', unit='inches')  # width of a capital M in inches
xlabheight <- max(strwidth(xlabs, unit='inches'))/emw
print(xlabheight)
# 2em margins top and bottom, 1.5 inches per panel + half the panel size for the column text label, xlabheight*em inches for the x labels
figheight <- nxgroups*1.5 + 1.5/2 + 4*em + xlabheight*emw



# Loop over devices to save both pdf and svgs
devs <- list(svglite, pdf)
outs <- c(out.svg, out.pdf)
for (i in 1:2) {
    cat('making', outs[i], '\n')
    devs[[i]](file=outs[i], width=figwidth, height=figheight)

    # from above: each element on the x-axis is given 1.5 ems of space. Give the
    # left-side group labels 6 ems of space. Very rough, probably won't always work.
    layout(matrix(1:((nxgroups+1)*(nygroups+2)), nrow=(nxgroups+1), byrow=TRUE),
        heights=c(2,rep(4,nxgroups)), widths=c(6,6,rep(1.5*nxlabs, nygroups)))
    par(oma=c(2+xlabheight,2,2,2))

    # just plot labels for y groups
    textpane('')
    textpane('')
    for (yg in ygroups) {
        textpane(strsplit(yg, '\n')[[1]])
    }

    # plots begin here
    for (xg in xgroups) {
        textpane(xg, legend.position='left')
        textpane('')

        # ranges based on error bars
        g <- d[x_group_factor_class == xg]
        print(c(g$enr.boot.0.95.lb, g$enr.boot.0.95.ub))
        # add in 0.95 and 1.05 to ensure we always have a reasonable amount
        # of spread around the 1 line.
        ylim <- range(c(0.95, 1.05, g$enr.boot.0.95.lb, g$enr.boot.0.95.ub), na.rm=TRUE)
        cat('uniform y-axis limit:', ylim, 'for group', yg, '\n')

        for (yg in ygroups) {
            par(mar=c(1/2,1/2,0,0))
            g <- g[y_group_factor_class == yg]
            meta.cols <- colnames(g)
            # for general BEDs (i.e., not qBEDs), this meta column is named 'quantile'
            # but the data are not necessarily quantiles. They could be data classes
            # with no numeric interpretation, like "inside a peak', a gene, or a
            # chromHMM state.
            meta.cols <- meta.cols[(1:length(meta.cols)) < which(meta.cols=='quantile')]
            signals <- g
print(signals)
            cat('plotting', nrow(signals), 'points\n')
            if (nrow(signals) == 0) {
                # Plot empty box placeholder
                plot(1, pch=NA, xaxt='n', yaxt='n')
            } else {
                color <- 'black'
                if (length(highlight) > 0)
                    color <- ifelse(signals[[highlight.key]] == highlight.value, 'green', 'black')
                plot.enrich(es=as.data.frame(signals), ylim=ylim, type='l', ltype='p',
                     yaxt=ifelse(yg==ygroups[1], 's', 'n'), lcol=color, xaxt='n')
                if (xg == tail(xgroups,1))
                    axis(side=1, at=1:nrow(signals), labels=signals[['x_label_factor_class']], las=3)
                abline(h=1, lty='dashed', col='grey')
            }
        }
    }
}

fwrite(d, file=out.csv)

if ('snakemake' %in% ls()) {
    sink()
}
