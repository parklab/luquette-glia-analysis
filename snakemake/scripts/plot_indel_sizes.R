#!/usr/bin/env Rscript

# detect script being run by snakemake
# if so, make a mock commandArgs function
if ('snakemake' %in% ls()) {
    logfile <- snakemake@log[[1]]
    con <- file(logfile, 'w')
    sink(con, type='output')
    sink(con, type='message')

    commandArgs <- function(...) unlist(c(
        snakemake@input[1:2], snakemake@output[1:3]
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 4) {
    cat('assumes in.rda contains only one object\n')
    stop("usage: plot_indel_sizes.R neuron_indels.rda oligo_indels.rda out.svg out.pdf out.csv")
}

nfile <- args[1]
gfile <- args[2]
outsvg <- args[3]
outpdf <- args[4]
outcsv <- args[5]

if (file.exists(outsvg))
    stop(paste('output file', outsvg, 'already exists, please delete it first'))
if (file.exists(outpdf))
    stop(paste('output file', outpdf, 'already exists, please delete it first'))
if (file.exists(outcsv))
    stop(paste('output file', outcsv, 'already exists, please delete it first'))

suppressMessages(library(data.table))
suppressMessages(library(extrafont))
suppressMessages(library(svglite))

if (!("Arial" %in% fonts()))
    stop("Arial font not detected; did you load extrafonts and run font_import() with the appropriate path?")


n <- get(load(nfile))
g <- get(load(gfile))
nsize <- nchar(n$altnt) - nchar(n$refnt)
gsize <- nchar(g$altnt) - nchar(g$refnt)

fwrite(table(rbind(data.frame(CellType='neuron', IndelLength=nsize), data.frame(CellType='oligo', IndelLength=gsize))), file=outcsv)

axisat <- c(min(c(nsize,gsize)), -20, -10, -5, -1, 1, 5, 10, max(c(nsize,gsize)))

figwidth=4
figheight=6
devs <- list(svglite, pdf)
outs <- c(outsvg, outpdf)
for (i in 1:2) {
    devs[[i]](file=outs[i], width=figwidth, height=figheight)

    xlim <- range(c(nsize, gsize))
    layout(1:2)
    plot(table(nsize), col=1, xlim=xlim, xaxt='n',
        main='Neurons', xlab='Indel size (< 0 = deletion)', ylab='Number')
    axis(side=1, at=axisat)
    plot(table(gsize), col=2, xlim=xlim, xaxt='n',
        main='Oligo', xlab='Indel size (< 0 = deletion)', ylab='Number')
    axis(side=1, at=axisat)
}

if ('snakemake' %in% ls()) {
    sink()
}
