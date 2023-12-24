#!/usr/bin/env Rscript

# detect script being run by snakemake
# if so, make a mock commandArgs function
if ('snakemake' %in% ls()) {
    logfile <- snakemake@log[[1]]
    con <- file(logfile, 'w')
    sink(con, type='output')
    sink(con, type='message')

    commandArgs <- function(...) unlist(c(
        snakemake@input['rda'],
        #snakemake@output['rda'],
        snakemake@output['csv'],
        snakemake@params['fdr']
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}


args <- commandArgs(trailingOnly=T)

if (length(args) != 3)
    stop(sprintf("usage: change_scan2_fdr.R scan2_input.rda out.csv new_fdr"))

infile <- args[1]
#out.rda <- args[2]
out.csv <- args[2]
fdr <- as.numeric(args[3])

if (fdr < 0 | fdr > 1)
    stop('0 < fdr <= 1 required')

#for (f in c(out.rda, out.csv)) {
for (f in c(out.csv)) {
    if (file.exists(f))
        stop(sprintf("output file %s already exists, please delete it first", f))
}


# DO WORK
suppressMessages(library(scan2))

load(infile)  # loads 'results'

# not generally OK, but we don't use this in the tiny use-case for changed-FDR objects
results@spatial.sensitivity <- NULL   

results@gatk$sample <- results@single.cell
results@gatk$bulk.sample <- results@bulk
# can't have sample names in columns - rbind would like uniform column names
colnames(results@gatk)[colnames(results@gatk) == results@single.cell] <- 'scgt'

old.results <- results
results <- copy(old.results)   # data.table is weird

results <- call.mutations(results, target.fdr=fdr)

fwrite(results@gatk[pass == TRUE], file=out.csv)

#save(results, file=out.rda, compress=FALSE)

if ('snakemake' %in% ls()) {
    sink()
}
