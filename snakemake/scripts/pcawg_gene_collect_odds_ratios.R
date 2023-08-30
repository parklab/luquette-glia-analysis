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
        snakemake@input['cancer_ors']
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 2) {
    stop('usage: pcawg_gene_collect_odds_ratios.R out.csv cancer_odds_ratios1.rda [ cancer_odd_ratios2.rda ... cancer_odds_ratiosN.rda ]')
}


out.csv <- args[1]
cancer.fs <- args[-(1)]

for (f in c(out.csv)) {
    if (file.exists(f))
        stop(paste('output file', f, 'already exists, please delete it first'))
}

suppressMessages(library(data.table))

data <- lapply(cancer.fs, function(f) {
    load(f)  # among other things, loads 'tops', 'types' and 'fts'
    list(tops=tops, types=types, odds=fts[1,])
})
odds.mat <- sapply(1:length(data), function(i) data[[i]]$odds)
colnames(odds.mat) <- sapply(data, function(d) d$types)
types <- colnames(odds.mat)
xaxis <- data[[1]]$tops

odds.mat <- data.table(TopNGenes=xaxis, odds.mat)[order(TopNGenes)]
fwrite(odds.mat, file=out.csv)

if ('snakemake' %in% ls()) {
    sink()
}
