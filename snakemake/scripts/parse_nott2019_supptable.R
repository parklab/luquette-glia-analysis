#!/usr/bin/env Rscript

# detect script being run by snakemake
# if so, make a mock commandArgs function
if ('snakemake' %in% ls()) {
    logfile <- snakemake@log[[1]]
    con <- file(logfile, 'w')
    sink(con, type='output')
    sink(con, type='message')

    commandArgs <- function(...) unlist(c(
        snakemake@input['xlsx'],
        snakemake@params['out_prefix']
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 2) {
    stop('usage: parse_nott2019_supptable.R supplementary_table5.xlsx out_prefix')
}

in.xlsx <- args[1]
out.prefix <- args[2]

celltypes <- c('astrocyte', 'microglia', 'neuron', 'oligo')
rtypes <- c('enhancer', 'promoter')

# IMPORTANT: worksheet number -> type
worksheet.number <- c(
    astrocyte_enhancer=5, astrocyte_promoter=6,
    neuron_enhancer=7, neuron_promoter=8,
    oligo_enhancer=9, oligo_promoter=10,
    microglia_enhancer=11, microglia_promoter=12
)

# this must match naming scheme above for worksheet.number
file.roots <- as.character(outer(celltypes, rtypes, paste, sep='_'))

tmp.out.files <- paste0(out.prefix, '/', file.roots, '.csv')
out.files <- paste0(out.prefix, '/', file.roots, '.bed')
for (f in c(out.files, tmp.out.files))
    if (file.exists(f))
        stop(paste0('output file ', f, ' already exists, please delete it first'))

suppressMessages(library(scan2))
suppressMessages(library(rio))


for (i in 1:length(file.roots)) {
    tmp.csv <- tmp.out.files[i]
    out.csv <- out.files[i]
    file.root <- file.roots[i]
    rio::convert(in_file=in.xlsx, out_file=tmp.csv, in_opts=list(sheet=worksheet.number[file.root]))
    # first two lines are title and blank
    tmp <- data.table::fread(tmp.csv, skip=2, col.names=c('chr','start','end'))
    # file contains duplicated lines
    tmp <- tmp[!duplicated(paste(chr,start,end))]
    # must be column 4
    tmp$type <- file.root
    fwrite(tmp, file=out.csv, sep='\t', col.names=FALSE)
}

if ('snakemake' %in% ls()) {
    sink()
}
