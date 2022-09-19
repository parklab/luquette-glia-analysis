#!/usr/bin/env Rscript

# detect script being run by snakemake
# if so, make a mock commandArgs function
if ('snakemake' %in% ls()) {
    logfile <- snakemake@log[[1]]
    con <- file(logfile, 'w')
    sink(con, type='output')
    sink(con, type='message')

    commandArgs <- function(...) unlist(c(
        snakemake@input[c('neuron','oligo')], snakemake@output['csv']
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 3) {
    stop("usage: rda_to_csv.r neuron_scores.csv oligo_scores.csv out.csv")
}

nscore.csv <- args[1]
gscore.csv <- args[2]
out.csv <- args[3]

if (file.exists(out.csv))
    stop(paste('outut file', out.csv, 'already exists, please delete it first'))

suppressMessages(library(scan2))

n <- fread(nscore.csv)
g <- fread(gscore.csv)

n$AgeSignifAdj <- -log10(p.adjust(10^-n$AgeSignif))
n[, SigIncluded := (AgeSignifAdj >= 2 & PercentContribution >= 1) | PercentStepwiseMeanResidualReduction > 50]

g$AgeSignifAdj <- -log10(p.adjust(10^-g$AgeSignif))
g[, SigIncluded := (AgeSignifAdj >= 2 & PercentContribution >= 1) | PercentStepwiseMeanResidualReduction > 50]

combined <- rbind(n,g)

fwrite(combined, file=out.csv)

if ('snakemake' %in% ls()) {
    sink()
}
