#!/usr/bin/env Rscript

# detect script being run by snakemake
# if so, make a mock commandArgs function
if ('snakemake' %in% ls()) {
    logfile <- snakemake@log[[1]]
    con <- file(logfile, 'w')
    sink(con, type='output')
    sink(con, type='message')

    commandArgs <- function(...) unlist(c(
        snakemake@input[1], snakemake@output[1]
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 2) {
    stop("usage: make_cosmic_indel_table.R in_cosmic.csv out.csv")
}


cosmic.csv <- args[1]
outcsv <- args[2]

if (file.exists(outcsv))
    stop(paste('output file', outcsv, 'already exists, please delete it first'))

library(scan2)

# N.B.
# Unlike COSMIC SBS, the indel database has no " -    " values
cosmic <- fread(cosmic.csv, header=T, stringsAsFactors=F)
head(cosmic)


id83 <- levels(id83(c()))

# In the COSMIC database, ID83 channel names are formatted differently
# from SigProfilerMatrixGenerator (spmgr). E.g.:
# COSMIC
# [1] "DEL_C_1_0"  "DEL_C_1_1"  "DEL_C_1_2"  "DEL_C_1_3"  "DEL_C_1_4" 
# [6] "DEL_C_1_5+" ...
# SPMGR
# [1] "1:Del:C:0" "1:Del:C:1" "1:Del:C:2" "1:Del:C:3" "1:Del:C:4" "1:Del:C:5"
# This function translates COSMIC names -> SPMGR names
id83.cosmic.to.spmgr <- function(x) {
    sapply(strsplit(sub('MH', 'M', gsub('\\+', '', gsub('_', ':', sub('repeats', 'R', sub('DEL', 'Del', sub('INS', 'Ins', x)))))), ':'), function(elts) paste(elts[3], elts[1], elts[2], elts[4], sep=':'))
}

cosmic[[1]] <- id83.cosmic.to.spmgr(cosmic[[1]])
colnames(cosmic)[1] <- 'MutType'
setkey(cosmic, MutType)
cosmic <- cosmic[id83]
head(cosmic)

fwrite(cosmic, file=outcsv)

if ('snakemake' %in% ls()) {
    sink()
}
