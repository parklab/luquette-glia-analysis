#!/usr/bin/env Rscript

# detect script being run by snakemake
# if so, make a mock commandArgs function
if ('snakemake' %in% ls()) {
    #logfile <- snakemake@log[[1]]
    #con <- file(logfile, 'w')
    #sink(con, type='output')
    #sink(con, type='message')

    commandArgs <- function(...) unlist(c(
        snakemake@input[1], # rda digest file
        snakemake@output[1:2],  # svg/pdf output file names
        snakemake@params[1] # binsize: 100, 1k, 10k, 100k, 1m
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 4) {
    stop("usage: plot_chrom_depth.R digest.rda out.svg out.pdf binsize")
}

in.rda <- args[1]
out.svg <- args[2]
out.pdf <- args[3]
binsize <- args[4]


if (file.exists(out.svg))
    stop(paste('output file', out.svg, 'already exists, please delete it first'))
if (file.exists(out.pdf))
    stop(paste('output file', out.pdf, 'already exists, please delete it first'))

suppressMessages(library(extrafont))
suppressMessages(library(svglite))
suppressMessages(library(GenomicRanges))
#suppressMessages(library(BSgenome))
#suppressMessages(library(BSgenome.Hsapiens.UCSC.hg19))

if (!("Arial" %in% fonts()))
    stop("Arial font not detected; did you load extrafonts and run font_import() with the appropriate path?")

load(in.rda)

figwidth <- 9 # inches
figheight <- 4 
# Loop over devices to save both pdf and svgs
devs <- list(svglite, pdf)
outs <- c(out.svg, out.pdf)
for (i in 1:2) {
    devs[[i]](file=outs[i], width=figwidth, height=figheight)
    par(mar=c(5,4,1,1))
    chrom <- seqnames(gs[[1]])[1]
    pos <- apply(ps[[binsize]], 1, mean) # get the bin midpoint
    readmat <- ys[[binsize]]
    matplot(pos, 100*t(t(readmat) / colSums(readmat)),
        type='l', lty='solid', col=rgb(0,0,0,alpha=0.1), xlim=c(0, max(pos)),
        ylab='Percent reads per bin', xlab=paste(chrom, 'position'))
}
