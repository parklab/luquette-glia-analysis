#!/usr/bin/env Rscript

# detect script being run by snakemake
# if so, make a mock commandArgs function
if ('snakemake' %in% ls()) {
    logfile <- snakemake@log[[1]]
    con <- file(logfile, 'w')
    sink(con, type='output')
    sink(con, type='message')

    commandArgs <- function(...) unlist(c(
        snakemake@params['sigtype'],
        snakemake@params['plot_title'],
        snakemake@input['csv'],
        snakemake@output['svg'],
        snakemake@output['pdf'],
        snakemake@output['csv']
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 6) {
    stop("usage: plot_raw_spectrum.R {sbs96,id83} plot_title mutations.csv out.svg out.pdf out.csv")
}

sigtype <- args[1]
plot.title <- args[2]
muts.csv <- args[3]
outsvg <- args[4]
outpdf <- args[5]
outcsv <- args[6]

if (!(sigtype %in% c('sbs96','id83')))
    stop(paste('sigtype must be one of "sbs96" or "id83", case sensitive. got', sigtype))

if (file.exists(outsvg))
    stop(paste('output file', outsvg, 'already exists, please delete it first'))
if (file.exists(outpdf))
    stop(paste('output file', outpdf, 'already exists, please delete it first'))
if (file.exists(outcsv))
    stop(paste('output file', outcsv, 'already exists, please delete it first'))

suppressMessages(library(scan2))
suppressMessages(library(data.table))
suppressMessages(library(svglite))

muts <- fread(muts.csv)
plotf <- NULL
if (sigtype == 'sbs96') {
    muts$mutsig <- sbs96(muts$mutsig)
    plotf <- plot.sbs96
} else if (sigtype == 'id83') {
    muts$mutsig <- id83(muts$mutsig)
    plotf <- plot.id83
} else {
    stop(paste('unrecognized sigtype =', sigtype))
}


figwidth=7
figheight=2
devs <- list(svglite, pdf)
outs <- c(outsvg, outpdf)
for (i in 1:2) {
    devs[[i]](file=outs[i], width=figwidth, height=figheight)
    par(mar=c(1,4,3,1))
    s <- as.spectrum(muts$mutsig, eps=0, fraction=TRUE)
    plotf(x=1, spectrum=s, main=plot.title)
}

d <- data.table(MutType=names(s), Spectrum=as.vector(s))
fwrite(d, file=outcsv)

if ('snakemake' %in% ls()) {
    sink()
}
