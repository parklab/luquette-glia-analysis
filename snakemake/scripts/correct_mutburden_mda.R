#!/usr/bin/env Rscript

# detect script being run by snakemake
# if so, make a mock commandArgs function
if ('snakemake' %in% ls()) {
    logfile <- snakemake@log[[1]]
    con <- file(logfile, 'w')
    sink(con, type='output')
    sink(con, type='message')

    if (!('prevburdens' %in% names(snakemake@input)))
        snakemake@input['prevburdens'] <- 'Notafile'

    commandArgs <- function(...) unlist(c(
        snakemake@input['mutburden'],
        snakemake@input['expomat'],
        snakemake@output['csv']
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 3) {
    cat("expomat.csv is the matrix of signature exposures. It *MUST* have been computed with Signature B present or else this correction will not succeed.\n")
    cat("Furthermore, signature B exposure should *NOT* be present in expomat.csv.\n")
    stop("usage: correct_mutburden_mda.R mutburden.csv expomat.csv out.csv")
}

mutburden.file <- args[1]
expomat.file <- args[2]
out.csv <- args[3]

if (file.exists(out.csv))
    stop(paste('output file', out.csv, 'already exists, please delete it first'))

suppressMessages(library(scan2))

mutburden <- fread(mutburden.file)
expomat <- fread(expomat.file)

expo.sums <- colSums(expomat[,-1])
expo.sums <-data.table(sample=names(expo.sums), genome.burden=expo.sums)

# Rename the current genome.burden column
mutburden[, uncorrected.genome.burden := genome.burden]
mutburden$genome.burden <- NULL

# Corrected genome burden now occupies the "genome.burden" column.
mutburden <- mutburden[expo.sums,,on=.(sample)]

fwrite(mutburden, file=out.csv)

if ('snakemake' %in% ls()) {
    sink()
}
