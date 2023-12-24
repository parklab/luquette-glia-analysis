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
        snakemake@input['meta'],
        snakemake@input[['csvs']]
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}


args <- commandArgs(trailingOnly=T)

if (length(args) < 3) {
    cat("shared_count.csv files only contain pair-pair comparisons within an individual\n")
    stop("usage: figX_number_snvs.R out.pdf metadata.csv shared_counts1.csv [ shared_counts2.csv ... shared_countsN.csv ]")
}

out.pdf <- args[1]
meta.file <- args[2]
count.files <- args[-(1:2)]

for (f in c(out.pdf)) {
    if (file.exists(f))
        stop(sprintf("output file %s already exists, please delete it first", f))
}

# DO WORK
library(data.table)

meta <- fread(meta.file)
x <- rbindlist(lapply(count.files, fread))

# joining on sample1 is only valid assuming pairs are only constructed within an individual
x <- meta[x,, on=.(sample=sample1)]
x <- meta[x,, on=.(sample=sample2)]
# column name collisions: sample1 column names are preceeded by "i."

# Get only OL-OL pairs
y <- x[type=='oligo' & i.type=='oligo', .(age, donor, sample, i.sample, n.shared)]

ages <- sapply(split(y$age, y$donor), head, 1)
ns <- split(y$n.shared, y$donor)[order(ages)]

pdf(file=out.pdf, width=2.5, height=1.75, pointsize=8)
par(mar=c(4,4,2,0.1))
stripchart(ns, vertical=T, pch=20, cex=1/2, method='jitter', jitter=1/4, las=3, ylab='Number of shared sSNVs', xlab='Subject ID')
dev.off()

if ('snakemake' %in% ls()) {
    sink()
}
