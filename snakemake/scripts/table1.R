#!/usr/bin/env Rscript

# detect script being run by snakemake
# if so, make a mock commandArgs function
if ('snakemake' %in% ls()) {
    logfile <- snakemake@log[[1]]
    con <- file(logfile, 'w')
    sink(con, type='output')
    sink(con, type='message')

    commandArgs <- function(...) unlist(c(
        snakemake@input['meta'],
        snakemake@input['metrics'],
        snakemake@output['cells'],
        snakemake@output['muts'],
        snakemake@params[['samples']]
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}


args <- commandArgs(trailingOnly=T)

if (length(args) < 5)
    stop("usage: table1.R metadata.csv metrics.csv cell_table.csv mutation_table.csv sampleID_1 [ sampleID_2 ... sampleID_N ]")

meta.file <- args[1]
metrics.file <- args[2]
cells.out.csv <- args[3]
muts.out.csv <- args[4]
samples <- args[-(1:4)]

for (f in c(cells.out.csv, muts.out.csv)) {
    if (file.exists(f))
        stop(sprintf("output file %s already exists, please delete it first", f))
}

suppressMessages(library(scan2))

meta <- fread(meta.file)[sample %in% samples]
meta[, ageclass := factor(ageclass, levels=c('infant','adolescent','adult','elderly'), ordered=T)]

metrics <- fread(metrics.file)[sample %in% samples]

t1a <- meta[, .(Age=age[1], `Neurons`=sum(type=='neuron'), `OLs`=sum(type=='oligo'), ageclass=ageclass[1]),
    by=donor][order(Age)]

fwrite(t1a, file=cells.out.csv)

t1b <- meta[metrics,,on=.(sample)][,
    # amp != 'MDA': SCAN2 rescue calls should be ignored for MDA. They were
    # produced due to uniform running of the pipeline.
    .(`Mutations`=sum(calls.vafonly+1*(amp != 'MDA')*calls.rescue)),
    by=.(type,ageclass,muttype)][order(ageclass)]

fwrite(t1b, file=muts.out.csv)

if ('snakemake' %in% ls()) {
    sink()
}
