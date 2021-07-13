#!/usr/bin/env Rscript

# detect script being run by snakemake
# if so, make a mock commandArgs function
if ('snakemake' %in% ls()) {
    logfile <- snakemake@log[[1]]
    con <- file(logfile, 'w')
    sink(con, type='output')
    sink(con, type='message')

    commandArgs <- function(...) unlist(c(
        snakemake@input[1],      # table from make_roadmap_enrichment_table.R
        snakemake@output[1:3]   # output svg/pdf/csv ||no longer: /jpeg
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 4) {
    stop("usage: make_roadmap_enrichment_grid.R enrichment_table.csv out.svg out.pdf out.csv")
}

in.csv <- args[1]
out.svg <- args[2]
out.pdf <- args[3]
out.csv <- args[4]

if (file.exists(out.svg))
    stop(paste('output file', out.svg, 'already exists, please delete it first'))
if (file.exists(out.pdf))
    stop(paste('output file', out.pdf, 'already exists, please delete it first'))
if (file.exists(out.csv))
    stop(paste('output file', out.csv, 'already exists, please delete it first'))

library(data.table)

d <- fread(in.csv)

options(scipen=1000)
binsizes <- as.character(sort(as.integer(unique(d$binsize))))
nbins <- length(binsizes)
print(binsizes)

marks <- sort(unique(d$mark))
nmarks <- length(marks)
print(marks)

eids <- sort(unique(d$eid))
brain.eids <- sort(unique(d$eid[d$GROUP == 'Brain']))
neids <- length(eids)
print(eids)
cat('brain only\n')
print(brain.eids)

# just make a blank panel with text 'txt' in the middle
textpane <- function(txt) {
    par(mar=c(0,0,0,0))
    plot(1, pch=NA, bty='n', xaxt='n', yaxt='n')
    legend('center', legend=txt, bty='n', cex=1.4)
}

suppressMessages(library(extrafont))
suppressMessages(library(svglite))


if (!("Arial" %in% fonts()))
    stop("Arial font not detected; did you load extrafonts and run font_import() with the appropriate path?")

# get height (in lines) of longest class name (will be in x-axis label)
em <- strheight('M', unit='inches')  # height of a capital M in inches
# 4 - y axis, 1 - right margin, 1.5 inches per panel (plus text label panel)
figwidth <- (4 + 1)*em + (1+nmarks)*1.5  # in inches
# 2em margins top and bottom, 1.5 inches per panel + 0.5 inch for the column text label
figheight <- nbins*1.5 + 0.5 + 2*em


# To start, just look at the prefrontal cortex data
#d <- d[eid=='E073',]
# is.na: NA quantiles are the excluded and outside regions
ylim <- range(d$enr[!is.na(d$quantile)], na.rm=TRUE)
cat('uniform y-axis limit:', ylim, '\n')

# Loop over devices to save both pdf and svgs
devs <- list(svglite, pdf) #, jpeg)
outs <- c(out.svg, out.pdf) #, out.jpeg)
for (i in 1:2) { #3) {
    cat('making', outs[i], '\n')
    devs[[i]](file=outs[i], width=figwidth, height=figheight)
    # jpeg() defaults to units in pixels
    #if (i == 3) {
        #devs[[i]](file=outs[i], width=figwidth, height=figheight, unit="in", res=600)
    #} else {
        #devs[[i]](file=outs[i], width=figwidth, height=figheight)
    #}

    layout(matrix(1:((nbins+1)*(nmarks+1)), nrow=(nbins+1)),
        heights=c(1,rep(4,nbins)))
    par(oma=c(2,2,2,2))
    for (bs in c('', binsizes)) {
        textpane(bs)
    }
    for (m in marks) {
        textpane(m)
        for (bs in binsizes) {
            par(mar=c(1/2,1/2,0,0))
            for (e in eids) {
                plotf <- if (e == eids[1]) plot else lines
                dd <- d[mark == m & binsize==bs & eid == e]
                plotf(dd$quantile, dd$enr, ylim=ylim,
                    type='l',
                    #type='b', pch=20,
                    xaxt=ifelse(bs==tail(binsizes,1), 's', 'n'),
                    yaxt=ifelse(m==marks[1], 's', 'n'))
            }
            for (e in brain.eids) {
                dd <- d[mark == m & binsize==bs & eid == e]
                lines(dd$quantile, dd$enr, col=3,
                    type='l')
                    #type='b', pch=20)
            }
            abline(h=1, lty='dashed', col='grey')
        }
    }
}

fwrite(d, file=out.csv)

if ('snakemake' %in% ls()) {
    sink()
}
