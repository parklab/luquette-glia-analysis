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
if (length(args) != 5) {
    cat('assumes in.rda contains only one object\n')
    stop("usage: plot_indel_sizes.R neuron_indels.csv oligo_indels.csv out.svg out.pdf out.csv")
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
suppressMessages(library(svglite))


n <- fread(nfile)
g <- fread(gfile)
nsize <- nchar(n$altnt) - nchar(n$refnt)
gsize <- nchar(g$altnt) - nchar(g$refnt)

fwrite(table(rbind(data.frame(CellType='neuron', IndelLength=nsize), data.frame(CellType='oligo', IndelLength=gsize))), file=outcsv)

axisat <- c(min(c(nsize,gsize)), -20, -10, -5, -1, 1, 5, 10, max(c(nsize,gsize)))

figwidth=2.5
figheight=2.5
devs <- list(svglite, pdf)
outs <- c(outsvg, outpdf)
for (i in 1:2) {
    devs[[i]](file=outs[i], width=figwidth, height=figheight, pointsize=5)

    #xlim <- range(c(gsize, nsize))
    xlim <- c(-30, 30)
    gdens <- density(gsize)
    ndens <- density(nsize)
    ylim <- range(pretty(range(c(0, ndens$y, gdens$y))))
    plot(ndens, col=1, ylim=ylim, xlim=xlim,
        main='Indel size distribution',
        xlab='Indel size (< 0 = deletion)', ylab='Density')
    lines(gdens, col=2)
    legend('topright', bty='n', legend=c('Neurons', 'Oligodendrocytes'), lwd=2, col=1:2)
}

if ('snakemake' %in% ls()) {
    sink()
}
