#!/usr/bin/env Rscript

# detect script being run by snakemake
# if so, make a mock commandArgs function
if ('snakemake' %in% ls()) {
    logfile <- snakemake@log[[1]]
    con <- file(logfile, 'w')
    sink(con, type='output')
    sink(con, type='message')

    commandArgs <- function(...) unlist(c(
        snakemake@output['pdf'],
        snakemake@input[['csvs']]
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}


args <- commandArgs(trailingOnly=T)

if (length(args) < 2) {
    cat("inX.csv files must be named {donor_ID}.csv!\n")
    stop("usage: fig3_maybe_number_shared_snvs.R out.pdf in1.csv [ in2.csv ... inN.csv]")
}

out.pdf <- args[1]
in.csvs <- args[-1]

for (f in c(out.pdf)) {
    if (file.exists(f))
        stop(sprintf("output file %s already exists, please delete it first", f))
}

# DO WORK
library(data.table)

x <- rbindlist(lapply(in.csvs, function(f) {
    print(f)
    x <- fread(f)
    cbind(donor=sub('.csv', '', sapply(strsplit(f, split='/'), tail, 1)), x)
}))
print(x)

pdf(file=out.pdf, height=2, width=1/2 + length(in.csvs)*strwidth('M', unit='inches'), pointsize=7)
par(mar=c(7,4,1,1))
stripchart(split(x$n.shared, x$donor),
    las=3, vertical=T, pch=20, method='jitter',
    cex=3/4, jitter=0.2, bty='l', ylab='Number of shared sSNVs')
dev.off()

if ('snakemake' %in% ls()) {
    sink()
}
