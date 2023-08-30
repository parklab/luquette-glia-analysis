#!/usr/bin/env Rscript

# detect script being run by snakemake
# if so, make a mock commandArgs function
if ('snakemake' %in% ls()) {
    logfile <- snakemake@log[[1]]
    con <- file(logfile, 'w')
    sink(con, type='output')
    sink(con, type='message')

print(snakemake@params['samples'])
    commandArgs <- function(...) unlist(c(
        snakemake@input['csv'],
        snakemake@params['qualtype'],
        snakemake@params['filter'],
        snakemake@output['csv'],
        unlist(snakemake@params['samples'])
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 4) {
    cat("This script is intended to subset the mutation tables output by scan2 rescue.\n")
    stop("usage: scan2_table_to_csv.R scan2_table.csv {A|AB|indel_A|indel_AB} {filtered|unfiltered} out.csv [ sample1 ... sampleN ]")
}

incsv <- args[1]
qualtype <- args[2]
filter <- tolower(args[3])
outcsv <- args[4]

samples <- NULL
if (length(args) > 4)
    samples <- args[-(1:4)]

if (!(qualtype %in% c('A', 'AB', 'indel_A', 'indel_AB')))
    stop('qualtype must be one of A, AB, indel_A or indel_AB (case sensitive)')

if (!(filter %in% c('filtered', 'unfiltered')))
    stop(paste('argument 3 must be either "filtered" or "unfiltered", got', filter))

if (file.exists(outcsv))
    stop(paste('output file', outcsv, 'already exists, please delete it first'))

suppressMessages(library(scan2))

muttab <- fread(incsv)

# outdated name scheme for mutation calls: A = VAF-based (pass), B = mutsig rescued (rescue)
# indel_ for indels; no prefix means SNV
mt <- ifelse(substr(qualtype, 1, 5) == 'indel', 'indel', 'snv')
add.rescue <- substr(qualtype, nchar(qualtype), nchar(qualtype)) == 'B'

muttab <- muttab[muttype == mt & (pass == TRUE | (add.rescue & rescue == TRUE))]

if (!is.null(samples))
    muttab <- muttab[sample %in% samples]

if (filter == 'filtered')
    muttab <- muttab[final.filter == FALSE]

fwrite(muttab, file=outcsv)

if ('snakemake' %in% ls()) {
    sink()
}
