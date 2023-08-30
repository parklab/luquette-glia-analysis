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
        snakemake@input['metadata'],
        snakemake@output['txt'],
        snakemake@params['min_reads_per_bulk']
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}


args <- commandArgs(trailingOnly=T)

if (length(args) != 4)
    stop(sprintf("usage: get_bulk_alignable_bps.R chr.rda metadata.csv out.txt min.reads.per.bulk"))

in.rda <- args[1]
meta.file <- args[2]
out.txt <- args[3]
min.reads.per.bulk <- as.integer(args[4])

for (f in c(out.txt)) {
    if (file.exists(f))
        stop(sprintf("output file %s already exists, please delete it first", f))
}

suppressMessages(library(scan2))

load(in.rda)  # loads 'tiles' among other things

meta <- fread(meta.file)

bulk.cols <- meta[amp == 'bulk']$sample

n.bulks <- length(bulk.cols)

n.basepairs <- sum(width(tiles[rowSums(mean.mat[,..bulk.cols] >= min.reads.per.bulk) == n.bulks,]))

fwrite(data.table(chr=as.character(seqnames(tiles)[1]),
    n.bulks=n.bulks,
    min.reads.per.bulk=min.reads.per.bulk,
    n.basepairs=n.basepairs), file=out.txt)

if ('snakemake' %in% ls()) {
    sink()
}
