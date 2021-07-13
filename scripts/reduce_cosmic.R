#!/usr/bin/env Rscript

# detect script being run by snakemake
# if so, make a mock commandArgs function
if ('snakemake' %in% ls()) {
    logfile <- snakemake@log[[1]]
    con <- file(logfile, 'w')
    sink(con, type='output')
    sink(con, type='message')

    commandArgs <- function(...) unlist(c(
        snakemake@input[1:2], snakemake@output[1]
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 3) {
    stop("usage: reduce_cosmic.R cosmic_scores.csv cosmic.csv out.csv")
}

scores.csv <- args[1]
cosmic.csv <- args[2]
out.csv <- args[3]

if (file.exists(out.csv))
    stop(paste('output file', out.csv, 'already exists, please delete it first'))

library(data.table)

scores <- fread(scores.csv)
cosmic <- fread(cosmic.csv)
head(cosmic)

sigs.to.keep <- scores[SigIncluded == TRUE, SigName]
cat("keeping", length(sigs.to.keep), "signatures:\n")
print(sigs.to.keep)
# column 1 is the mutation type column
cosmic <- cbind(cosmic[,1], cosmic[, ..sigs.to.keep])
head(cosmic)

fwrite(cosmic, file=out.csv)

if ('snakemake' %in% ls()) {
    sink()
}
