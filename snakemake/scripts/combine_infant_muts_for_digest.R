#!/usr/bin/env Rscript

# detect script being run by snakemake
# if so, make a mock commandArgs function
if ('snakemake' %in% ls()) {
    logfile <- snakemake@log[[1]]
    con <- file(logfile, 'w')
    sink(con, type='output')
    sink(con, type='message')

    commandArgs <- function(...) unlist(c(
        snakemake@output['csv'],
        snakemake@input[['csvs']]
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}


args <- commandArgs(trailingOnly=T)

if (length(args) < 2)
    stop(sprintf("usage: combine_infant_muts_for_digest.R out.csv in1.csv [ in2.csv ... inN.csv ]"))

outfile <- args[1]
infiles <- args[-1]

for (f in c(outfile)) {
    if (file.exists(f))
        stop(sprintf("output file %s already exists, please delete it first", f))
}

# DO WORK
library(data.table)

fwrite(rbindlist(lapply(infiles, fread)), file=outfile)

if ('snakemake' %in% ls()) {
    sink()
}
