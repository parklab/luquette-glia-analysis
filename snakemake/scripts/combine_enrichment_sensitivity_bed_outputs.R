#!/usr/bin/env Rscript

# detect script being run by snakemake
# if so, make a mock commandArgs function
if ('snakemake' %in% ls()) {
    logfile <- snakemake@log[[1]]
    con <- file(logfile, 'w')
    sink(con, type='output')
    sink(con, type='message')

    commandArgs <- function(...) unlist(c(
        snakemake@input['sens'],
        snakemake@input['depth'],
        snakemake@output['txt']
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}


args <- commandArgs(trailingOnly=T)

if (length(args) != 3)
    stop("usage: combine_enrichment_sensitivity_bed_outputs.R sensitivity.summary.txt depth.summary.txt out.txt")

sens.file <- args[1]
depth.file <- args[2]
outfile <- args[3]

for (f in c(outfile)) {
    if (file.exists(f))
        stop(sprintf("output file %s already exists, please delete it first", f))
}

suppressMessages(library(data.table))

# DO WORK
sens <- fread(sens.file)
depth <- fread(depth.file)
depth$width <- NULL   # sens and depth don't always agree on feature width

fwrite(sens[depth,,on=.(metaline,feature,sample)], file=outfile)

if ('snakemake' %in% ls()) {
    sink()
}
