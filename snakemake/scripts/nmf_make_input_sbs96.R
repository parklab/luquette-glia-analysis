#!/usr/bin/env Rscript

# detect script being run by snakemake
# if so, make a mock commandArgs function
if ('snakemake' %in% ls()) {
    logfile <- snakemake@log[[1]]
    con <- file(logfile, 'w')
    sink(con, type='output')
    sink(con, type='message')

    commandArgs <- function(...) unlist(c(
        snakemake@input['mat'], snakemake@output['converted_mat']
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}


args <- commandArgs(trailingOnly=T)

if (length(args) != 2)
    stop(sprintf("usage: make_input_sbs96.r input.all output.csv"))

infile <- args[1]
outfile <- args[2]

if (file.exists(outfile))
    stop(sprintf("output file %s already exists, please delete it first", outfile))


x <- read.table(infile,stringsAsFactors=F,header=T,check.names=F)

y <- data.frame(
    MutationType=substr(x$MutationType, 3,5),
    Trinucleotide=paste0(substr(x$MutationType,1,1),
                         substr(x$MutationType,3,3),
                         substr(x$MutationType,7,7)),
    x[,-1],
    check.names=F, stringsAsFactors=F)

colnames(y)[-(1:2)] <- paste0('NormalTissue::', colnames(y)[-(1:2)])
write.csv(y, file=outfile, row.names=FALSE)

if ('snakemake' %in% ls()) {
    sink()
}
