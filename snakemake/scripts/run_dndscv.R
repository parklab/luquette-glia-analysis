#!/usr/bin/env Rscript

# detect script being run by snakemake
# if so, make a mock commandArgs function
if ('snakemake' %in% ls()) {
    logfile <- snakemake@log[[1]]
    con <- file(logfile, 'w')
    sink(con, type='output')
    sink(con, type='message')

    commandArgs <- function(...) unlist(c(
        snakemake@input['muts'],
        snakemake@params['sampletag'],
        snakemake@output['rda'],
        snakemake@output['csv']
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 4) {
    cat('sampletag will be saved in the .rda output so that later scripts can propagate the tumor type or celltype_qualtype\n')
    stop("usage: run_dndscv.R muts.csv sampletag out.rda out.csv")
}

muts.csv <- args[1]
sampletag <- args[2]
outrda <- args[3]
outcsv <- args[4]

if (file.exists(outrda))
    stop(paste('output file', outrda, 'already exists, please delete it first'))
if (file.exists(outcsv))
    stop(paste('output file', outcsv, 'already exists, please delete it first'))

suppressMessages(library(data.table))
suppressMessages(library(dndscv))

muts <- fread(muts.csv)
# luckily these column names are the same for our somatic mut tables and cancer snvs
setcolorder(muts, c('sample', 'chr', 'pos', 'refnt', 'altnt'))
colnames(muts)[1:5] <- c('sampleID', 'chr', 'pos', 'ref', 'mut')  # colnames expected by dndscv

print(head(muts))

dnds <- dndscv(muts[,1:5])

save(dnds, sampletag, file=outrda)
fwrite(dnds$sel_cv, file=outcsv)

if ('snakemake' %in% ls()) {
    sink()
}
