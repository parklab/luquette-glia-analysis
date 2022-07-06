#!/usr/bin/env Rscript

# detect script being run by snakemake
# if so, make a mock commandArgs function
if ('snakemake' %in% ls()) {
    logfile <- snakemake@log[[1]]
    con <- file(logfile, 'w')
    sink(con, type='output')
    sink(con, type='message')

    commandArgs <- function(...) unlist(c(
        snakemake@input[1:4], snakemake@output[1:3]
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 7)
    stop("usage: barplot_enrich_2x2 nmut.enrich.rda ni.enrich.rda gmut.enrich.rda gi.enrich.rda output.svg output.pdf output.tsv")

nefile <- args[1]
neifile <- args[2]
gefile <- args[3]
geifile <- args[4]
outsvg <- args[5]
outpdf <- args[6]
outtsv <- args[7]


if (file.exists(outsvg))
    stop(paste('outut file', outsvg, 'already exists, please delete it first'))
if (file.exists(outpdf))
    stop(paste('outut file', outpdf, 'already exists, please delete it first'))
if (file.exists(outtsv))
    stop(paste('outut file', outtsv, 'already exists, please delete it first'))

suppressMessages(library(GenomicRanges))
suppressMessages(library(extrafont))
suppressMessages(library(svglite))
suppressMessages(library(mutenrich))


if (!("Arial" %in% fonts()))
    stop("Arial font not detected; did you load extrafonts and run font_import() with the appropriate path?")

ne <- read.e(nefile)
nei <- read.e(neifile)
ge <- read.e(gefile)
gei <- read.e(geifile)


# ISSUES WITH GENERAL APPLICABILITY:
#   - need to specify par(mar), usually for x labels of classes
#   - need to specify class collapsing
#   - need to specify plot sizes

# get height (in lines) of longest class name (will be in x-axis label)
em <- strheight('M', unit='inches')  # height of a capital M in inches
xheight <- max(strwidth(names(ne$real.obs), unit='inches'))/em  # size in lines, for par(mar)

# 4 - y axis, 1 - right margin
figwidth <- 2*((4 + 1)*em + length(ne$real.obs)*em*2)  # in inches
# give 3.5 inches for the plot area, xheight for x-axis and 2 lines for title
# xheight*em = size in INCHES instead of lines
figheight <- 2*(3.5 + xheight*em + 2*em)

# Loop over devices to save both pdf and svgs
devs <- list(svglite, pdf)
outs <- c(outsvg, outpdf)
for (i in 1:2) {
devs[[i]](file=outs[i], width=figwidth, height=figheight)
# use x11 for some debugging
#x11(width=figwidth, height=figheight)
#svglite(file=outsvg)

layout(matrix(1:4, nrow=2, byrow=T))
par(mar=c(xheight, 4, 2, 1))

dlist <- list(`Neuron SNVs`=ne, `Oligo SNVs`=ge, `Neuron indels`=nei, `Oligo indels`=gei)

# row1: snvs
ylim <- range(c(get.ylim(ne, shift=1, scale=100, use.boot=T), get.ylim(ge, shift=1, scale=100, use.boot=T)))
barplot.enrich(dlist[[1]], main=names(dlist)[1], barcol='black', ylim=ylim, family='Arial')
barplot.enrich(dlist[[2]], main=names(dlist)[2], barcol='red', ylim=ylim, family='Arial')

# row2: indels
ylim <- range(c(get.ylim(nei, shift=1, scale=100, use.boot=T), get.ylim(gei, shift=1, scale=100, use.boot=T)))
barplot.enrich(dlist[[3]], main=names(dlist)[3], barcol='black', ylim=ylim, family='Arial')
barplot.enrich(dlist[[4]], main=names(dlist)[4], barcol='red', ylim=ylim, family='Arial')

dev.off()
}


# Write table of statistics
stat.table <- do.call(rbind, lapply(1:length(dlist), function(i) {
    es <- esummary(dlist[[i]])
    label <- names(dlist)[i]
    cbind(label, class=rownames(es), es)
}))
write.table(stat.table, quote=F, sep='\t', row.names=F, col.names=T, file=outtsv)

if ('snakemake' %in% ls()) {
    sink()
}
