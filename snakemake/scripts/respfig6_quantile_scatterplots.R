#!/usr/bin/env Rscript

# detect script being run by snakemake
# if so, make a mock commandArgs function
if ('snakemake' %in% ls()) {
    logfile <- snakemake@log[[1]]
    con <- file(logfile, 'w')
    sink(con, type='output')
    sink(con, type='message')

    commandArgs <- function(...) unlist(c(
        snakemake@input['a'],
        snakemake@input['ab'],
        snakemake@output['pdf']
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}


args <- commandArgs(trailingOnly=T)

if (length(args) != 3)
    stop("usage: respfig6_quantile_scatter_plots.R consolidated_signals_A.csv consolidated_signals_AB.csv out.pdf")

# "A" and "AB" are the old ways of referring to VAF-based and VAF+rescue-based SCAN2 calls
a.file <- args[1]
ab.file <- args[2]
out.pdf <- args[3]

for (f in c(out.pdf)) {
    if (file.exists(f))
        stop(sprintf("output file %s already exists, please delete it first", f))
}

# DO WORK
library(data.table)

a <- fread(a.file)
ab <- fread(ab.file)

# Some sanity checks to make sure the tables are identical
if (!all(a$datasource == ab$datasource))
    stop("all datasource values must be the same")
if (!all(a$lineclass == ab$lineclass))
    stop("all lineclass values must be the same")
if (!all(a$muttype == ab$muttype))
    stop("all muttype values must be the same")
if (!all(a$mutsfrom == ab$mutsfrom))
    stop("all mutsfrom values must be the same")

pdf(file=out.pdf, width=2.5, height=2*2.5, pointsize=8)
layout(1:2)
for (mt in c('snv', 'indel')) {
    par(mar=c(4,4,1,1))
    plot(x=a[muttype==mt]$enr, y=ab[muttype==mt]$enr, pch=20, cex=1/2, log='xy',
        col=ifelse(a[muttype==mt]$mutsfrom=='pta_neuron', 'black', 'red'),
        xlab='Enrichment values excluding rescued calls',
        ylab='Enrichment values including rescued calls',
        bty='l', main=ifelse(mt=='snv', 'sSNVs', 'Indels'))
    abline(coef=0:1)
}
dev.off()

if ('snakemake' %in% ls()) {
    sink()
}
