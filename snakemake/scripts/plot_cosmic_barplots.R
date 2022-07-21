#!/usr/bin/env Rscript

# detect script being run by snakemake
# if so, make a mock commandArgs function
if ('snakemake' %in% ls()) {
    logfile <- snakemake@log[[1]]
    con <- file(logfile, 'w')
    sink(con, type='output')
    sink(con, type='message')

    # tags come from params
    commandArgs <- function(...) unlist(c(
        snakemake@input[1:2], snakemake@output[1:3]
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 5) {
    stop("usage: plot_cosmic_barplots.R expomat.csv mutburden.csv out.svg out.pdf out.csv")
}

expomat.csv <- args[1]
mutburden.csv <- args[2]
out.svg <- args[3]
out.pdf <- args[4]
out.csv <- args[5]

if (file.exists(out.svg))
    stop(paste('output file', out.svg, 'already exists, please delete it first'))
if (file.exists(out.pdf))
    stop(paste('output file', out.pdf, 'already exists, please delete it first'))

suppressMessages(library(data.table))
suppressMessages(library(extrafont))
suppressMessages(library(svglite))

if (!("Arial" %in% fonts()))
    stop("Arial font not detected; did you load extrafonts and run font_import() with the appropriate path?")


E <- fread(expomat.csv)
head(E)
m <- fread(mutburden.csv)  # needed for age/sorting

# Reorder E by increasing age
signames <- E$Sig
E <- E[,-1] # get rid of signature name column
colorder <- order(m$age)
E <- E[,..colorder]


# handle up to 12 signatures
sigcols <- c('#a6cee3','#1f78b4','#b2df8a',
    '#33a02c', '#fdbf6f','#ff7f00',
    '#fb9a99','#e31a1c',
    '#cab2d6','#6a3d9a','#ffff99', 'black')[1:nrow(E)]

# Sort the signatures from low to high contribution so that
# the largest contributors are the bars closest to the bottom
# of the plot.
# rearrange sigcols too so the same colors are used across plots
sigcols <- sigcols[order(rowSums(E), decreasing=T)]
signames <- signames[order(rowSums(E), decreasing=T)]
E <- E[order(rowSums(E), decreasing=T),]

# get height (in lines) of longest sample name (will be in x-axis label)
em <- strheight('M', unit='inches')  # height of a capital M in inches

xheight <- max(strwidth(colnames(E), unit='inches'))/em  # size in lines, for par(mar)
# give 3.5 inches for the plot area, xheight for x-axis and 2 lines for title
# xheight*em = size in INCHES instead of lines
figheight <- 2*(3.5 + xheight*em + 2*em)

# 4 - y axis, 1 - right margin
figwidth <- (4 + 1)*em + ncol(E)*em*1.2  # in inches

fwrite(t(E), file=out.csv, row.names=TRUE)

# Loop over devices to save both pdf and svgs
devs <- list(svglite, pdf)
outs <- c(out.svg, out.pdf)
for (i in 1:2) {
    devs[[i]](file=outs[i], width=figwidth, height=figheight)
    layout(1:2)
    par(mar=c(xheight,4,1,1))
    barplot(as.matrix(E), col=sigcols, las=3, border=NA, ylab='Number of mutations')
    legend('topleft', legend=signames, fill=sigcols)

    x <- apply(E, 2, function(col) col/sum(col))
    barplot(x, col=sigcols, las=3, border=NA, ylab='Fraction of mutations')
    dev.off()
}

if ('snakemake' %in% ls()) {
    sink()
}
