#!/usr/bin/env Rscript

# detect script being run by snakemake
# if so, make a mock commandArgs function
if ('snakemake' %in% ls()) {
    logfile <- snakemake@log[[1]]
    con <- file(logfile, 'w')
    sink(con, type='output')
    sink(con, type='message')

    commandArgs <- function(...) unlist(c(
        snakemake@input[1], snakemake@output[1]
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 2) {
    stop("usage: rda_to_csv.r in.csv out.csv")
}

incsv <- args[1]
outcsv <- args[2]

if (file.exists(outcsv))
    stop(paste('output file', outcsv, 'already exists, please delete it first'))

library(data.table)

x <- fread(incsv)

# order of columns is important: all metadata is expected to come
# BEFORE the quantile column, which signals the start of the enrichment
# summary object columns.
raw.quantile <- x[['quantile']]
idx <- which(colnames(x) == 'quantile')

# mutation signature class is considered a metacolumn so we can group
# by it in later plots.
metacols <- x[,1:(idx-1)]
metacols$sigclass <- sapply(strsplit(raw.quantile, "|||", fixed=TRUE), tail, 1)

enrobjcols <- x[,idx:ncol(x)]
enrobjcols$quantile <- sapply(strsplit(raw.quantile, "|||", fixed=TRUE), head, 1)

x <- cbind(metacols, enrobjcols)
fwrite(x, file=outcsv)

if ('snakemake' %in% ls()) {
    sink()
}
