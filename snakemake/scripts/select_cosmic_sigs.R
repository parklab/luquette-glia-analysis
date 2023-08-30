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
        snakemake@input
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 2) {
    stop("usage: select_cosmic_sigs.R out.csv scored_sigs1.csv [ scored_sigs2.csv ... scored_sigs3.csv ]")
}

out.csv <- args[1]
score.csvs <- args[-1]

if (file.exists(out.csv))
    stop(paste('outut file', out.csv, 'already exists, please delete it first'))

suppressMessages(library(scan2))

scores <- rbindlist(lapply(score.csvs, fread))

scores[, SigIncluded := (AgeSignifAdj >= 4 & PercentContribution >= 1) | PercentStepwiseMeanResidualReduction > 25]

fwrite(scores, file=out.csv)

if ('snakemake' %in% ls()) {
    sink()
}
