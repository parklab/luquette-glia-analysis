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
        snakemake@output['csv'],
        snakemake@input[['summaries']],
        snakemake@input[['skinny']]
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 4) {
    cat("WARNING: summaryI.mat files MUST correspond to skinnyI.mat files in the same order!\n")
    stop("usage: nmf_plot_selection_criteria.R output.pdf output.csv summary1.mat [ summary2.mat ... summaryN.mat ] skinny1.mat [ skinny2.mat ... skinnyN.mat ]")
}

outpdf <- args[1]
outcsv <- args[2]
remaining <- args[-(1:2)]
summary.files <- remaining[1:(length(remaining)/2)]
skinny.files <- remaining[(1+length(remaining)/2):length(remaining)]

if (file.exists(outpdf))
    stop(paste('output file', outpdf, 'already exists, please delete it first'))
if (file.exists(outcsv))
    stop(paste('output file', outcsv, 'already exists, please delete it first'))

suppressMessages(library(R.matlab))
suppressMessages(library(scan2))

data <- lapply(1:length(summary.files), function(i) {
    sumf <- summary.files[i]
    skif <- skinny.files[i]
    print(sumf)
    print(skif)
    m <- readMat(sumf)
    # Our NMF runs always check a single value of N for decomposing signatures.
    # Each run checks one N value, allowing multiple runs for parallelization.
    minsigs <- m$minSignatures[1,1]
    maxsigs <- m$maxSignatures[1,1]
    if (minsigs != maxsigs)
        stop(paste('expected only a single N value in matrix file', f))
    avgstab <- m$avgStability
    avgerrpct <- m$avgReconstructionErrorPercentage

    m <- readMat(skif)

    list(nsigs=minsigs, avgstab=avgstab, avgerrpct=avgerrpct, avgstab.per.proc=m$processStabAvg[1,])
})

vals <- data.table(NSignatures=sapply(data, function(x) x$nsigs),
                  AverageStability=sapply(data, function(x) x$avgstab),
                  AverageReconstructionErrorPct=sapply(data, function(x) x$avgerrpct))

# Add an extra matrix of per-signature average stabilities. Each row will have a different
# number of averages since each row indicates a different number of signatures.
min.n <- min(vals$NSignatures)
max.n <- max(vals$NSignatures)
stabs.per.proc <- sapply(1:length(data), function(i) {
    ret <- rep(NA, max.n)
    d <- data[[i]]
    ret[1:d$nsigs] <- d$avgstab.per.proc
    ret
})

vals <- cbind(vals, t(stabs.per.proc))

fwrite(vals, file=outcsv)

pdf(outpdf, height=9)
layout(matrix(1:2, ncol=1))
plot(vals$NSignatures, vals$AverageStability, type='b',
    xlab='Number of signatures', ylab='Average signature stability')
plot(vals$NSignatures, vals$AverageReconstructionErrorPct, type='b',
    xlab='Number of signatures', ylab='Average reconstruction error (%)')
dev.off()

if ('snakemake' %in% ls()) {
    sink()
}
