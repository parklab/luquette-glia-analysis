#!/usr/bin/env Rscript

# detect script being run by snakemake
# if so, make a mock commandArgs function
if ('snakemake' %in% ls()) {
    logfile <- snakemake@log[[1]]
    con <- file(logfile, 'w')
    sink(con, type='output')
    sink(con, type='message')

    commandArgs <- function(...) unlist(c(
        snakemake@input['mat'],
        snakemake@input['scan2_id83_correction'],
        snakemake@output['mat']
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}


args <- commandArgs(trailingOnly=T)
if (length(args) != 3)
    stop(sprintf("usage: make_input_id83.r input.all scan2_id83_correction.csv output.csv"))

infile <- args[1]
scan2.id83.correction.file <- args[2]
outfile <- args[3]

if (file.exists(outfile))
    stop(sprintf("output file %s already exists, please delete it first", outfile))

suppressMessages(library(scan2))

scan2.id83.correction <- fread(scan2.id83.correction.file)
setkey(scan2.id83.correction, muttype)

mat <- fread(infile, check.names=FALSE)

if (!all(sort(mat$MutationType) == sort(scan2.id83.correction$muttype))) {
    cat("mat$MutationType:\n")
    print(mat$MutationType)
    cat("scan2.id83.correction$muttype\n")
    print(scan2.id83.correction$muttype)
    stop('matrix and correction factor ID83 channel names do not match. Perhaps they use different nomenclature?')
}

# reorder (setkey)
scan2.id83.correction <- scan2.id83.correction[mat$MutationType]

for (i in 2:ncol(mat)) {
    total <- sum(mat[[i]], na.rm=TRUE)
    # if total was 0, just leave it alone (avoid div by 0)
    if (total > 0) {
        mat[[i]] <- mat[[i]] / scan2.id83.correction[[2]]
        # preserve total count and keep as integer
        mat[[i]] <- round(mat[[i]] * total/sum(mat[[i]], na.rm=TRUE))
    }
}

# SigProfiler matrices are tab separated
fwrite(mat, file=outfile, sep='\t')

if ('snakemake' %in% ls()) {
    sink()
}
