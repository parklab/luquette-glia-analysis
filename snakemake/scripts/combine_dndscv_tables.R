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
        paste0(snakemake@params[['sampletags']], '=', snakemake@input[['csvs']])
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 2) {
    cat('if "sampletag" is supplied, the column in the output will be named sampletag.\n')
    stop("usage: combine_dndscv_tables.R out.csv sampletag1=dnds_results1.csv [ sampletag2=dnds_results2.csv ... sampletagN=dnds_resultsN.csv ]")
}

outcsv <- args[1]

# parse the SAMPLETAG=FILE.csv argument format
spls <- strsplit(args[-1], '=')
dndscv.csvs <- sapply(spls, tail, 1)
# if strsplit returns only 1 string, then there was no '=' sign
names(dndscv.csvs) <- ifelse(sapply(spls, length) == 1, paste0('file_', 1:length(spls)), sapply(spls,head,1))

print(dndscv.csvs)

if (file.exists(outcsv))
    stop(paste('output file', outcsv, 'already exists, please delete it first'))

suppressMessages(library(data.table))
suppressMessages(library(dndscv))

# get a list of data.tables such that:
#  1. if the dnds results file is from a tumor analysis, the data.table
#     contains (gene, p-value of positive selection)
#  2. if the dnds results file is from a normal mutation analysis, the
#     data.table contains (gene, total exonic mutations)
# moreover, the first column will always be named gene_name and the
# second column will be the sampletag denoting which tumor/normal mutation
# dataset the value reflects.
dts <- lapply(1:length(dndscv.csvs), function(i) {
    csv <- dndscv.csvs[[i]]
    sampletag <- names(dndscv.csvs)[i]
    print(csv)
    #load(rda, verb=TRUE)
    print(sampletag)
    #dt <- as.data.table(dnds$sel_cv)
    dt <- fread(csv)
    setkey(dt, gene_name)
    if (substr(sampletag, 1, 7) == "cancer_") {
        # parsing cancer results - we only want the p-values so genes can be ranked
        # by selective signals.
        dt <- dt[,.(gene_name, pallsubs_cv)]
        # don't remove cancer_ prefix. it's useful to have in the table
        colnames(dt)[2] <- sampletag # substr(sampletag, 8, nchar(sampletag))
        return(dt)
    } else {
        # parsing normal mutation results - want the total number of mutations per
        # gene.
        dt[, n := n_syn + n_mis + n_non + n_spl]
        dt <- dt[,.(gene_name, n)]
        colnames(dt)[2] <- sampletag
        return(dt)
    }
})

newtab <- dts[[1]]
# just join all tables together on gene_name
for (i in 2:length(dts))
    newtab <- newtab[dts[[i]]]

fwrite(newtab, file=outcsv)

if ('snakemake' %in% ls()) {
    sink()
}
