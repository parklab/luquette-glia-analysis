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
    cat('assumes in.rda contains only one object\n')
    stop("usage: rda_to_csv.r in.rda out.csv")
}

inrda <- args[1]
outcsv <- args[2]

if (file.exists(outcsv))
    stop(paste('output file', outcsv, 'already exists, please delete it first'))

library(data.table)

x <- get(load(inrda))
fwrite(x, file=outcsv)

if ('snakemake' %in% ls()) {
    sink()
}
