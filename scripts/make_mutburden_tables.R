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
        snakemake@input['meta'],
        snakemake@params['muttype'],
        snakemake@input['prevburdens'],  # only for neurons; 'Notafile' for oligo
        snakemake@input['objects']       # variable length list
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 3) {
    stop("usage: make_mutburden_tables.R out.csv metadata.csv {snv|indel} {prevburdens.txt|Notafile} scan2_object1.rda [ scan2_object2.rda .. scan2_objectN.rda]")
}

out.csv <- args[1]
meta <- args[2]
muttype <- args[3]
prevburden <- args[4]
object.rdas <- args[-(1:4)]

if (file.exists(out.csv))
    stop(paste('output file', out.csv, 'already exists, please delete it first'))

if (muttype != 'snv' & muttype != 'indel')
    stop("muttype must be either 'snv' or 'indel', case sensitive")

suppressMessages(library(scan2))

meta <- fread(meta)[,.(sample,donor,age,plotage,color,type)]

objects <- lapply(object.rdas, function(rf) get(load(rf)))

dt <- data.table(sample=sapply(objects, function(o) o@single.cell),
                 nsom=sapply(objects, function(o) sum(o@mutburden[[muttype]]$ncalls)),
                 genome.burden=sapply(objects, function(o) o@mutburden$burden[2]))
if (prevburden != 'Notafile') {
    cat("Adding neuron burdens from previous paper..\n")
    dt <- rbind(dt, fread(prevburden)[,.(sample=sample,nsom=raw.calls, genome.burden=burden)])
}

meta <- meta[dt, on='sample']

fwrite(meta, file=out.csv)

if ('snakemake' %in% ls()) {
    sink()
}
