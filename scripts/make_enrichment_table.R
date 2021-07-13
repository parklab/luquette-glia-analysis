#!/usr/bin/env Rscript

# detect script being run by snakemake
# if so, make a mock commandArgs function
if ('snakemake' %in% ls()) {
    logfile <- snakemake@log[[1]]
    con <- file(logfile, 'w')
    sink(con, type='output')
    sink(con, type='message')

    commandArgs <- function(...) unlist(c(
        snakemake@output[1],     # output table
        snakemake@input[1],      # EID metadata from Roadmap 
        snakemake@input[2:length(snakemake@input)]  # variable sized list of esummary RDAs
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 2) {
    cat("WARNING: each esummary object must contain the same metadata\n")
    stop("usage: make_enrichment_table.R out.csv esummary1.rda [ esummary2.rda ... esummaryN.rda ]")
}

out.csv <- args[1]
esummary.rdas <- args[-1]

if (file.exists(out.csv))
    stop(paste('output file', out.csv, 'already exists, please delete it first'))

library(data.table)

full.table <- do.call(rbind, lapply(esummary.rdas, function(f) {
    load(f) # loads es and emeta
    do.call(rbind, lapply(names(es), function(ename) {
        m <- emeta[ename,] # meta columns to be duplicated for each row in es[[ename]]
        esum <- es[[ename]]
        # 'outside', 'excluded' and any other non-integer class will
        # generate warnings about NA, but these can be ignored. the
        # integers will be ordered properly followed by the non-integers
        # in their original order.
        esum <- esum[order(as.integer(rownames(esum))),]
        ret <- data.table(m, quantile=rownames(esum), esum)
       ret
    }))
}))

fwrite(full.table, file=out.csv)
