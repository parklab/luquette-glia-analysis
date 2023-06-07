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
    stop(sprintf("usage: make_input_id83.r input.all output.csv"))

infile <- args[1]
outfile <- args[2]

if (file.exists(outfile))
    stop(sprintf("output file %s already exists, please delete it first", outfile))


x <- read.table(infile,stringsAsFactors=F,header=T,check.names=F)

y <- data.frame(
    Type=substr(toupper(x$MutationType),3,5),
    Subtype=substr(x$MutationType,7,7),
    Indel_size=substr(x$MutationType,1,1),
    Repeat_MH_size=substr(x$MutationType,9,9),
    x[,-1],
    check.names=F, stringsAsFactors=F)

y$Subtype[y$Subtype=='R'] <- 'repeats'
y$Subtype[y$Subtype=='M'] <- 'MH'
y$Indel_size[y$Indel_size=='5'] <- '5+'
y$Repeat_MH_size[y$Repeat_MH_size=='5'] <- '5+'

colnames(y)[-(1:4)] <- paste0('NormalTissue::', colnames(y)[-(1:4)])
write.csv(y, file=outfile, row.names=FALSE)

if ('snakemake' %in% ls()) {
    sink()
}
