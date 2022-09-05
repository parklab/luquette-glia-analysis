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
        snakemake@params['sig'],
        snakemake@output['rda']
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 3) {
    cat('in.rda must contain only a single object. mutation .RData files and permutation .RData files are differentiated by checking this object for the data.table/data.frame or GenomicRanges/CompressedGRangesList classes, respectively. If the object is not one of these two classes, the script will exit.\n')
    cat('the object in in.rda must contain a column named `mutsig` that contains SBS96 channel info\n')
    cat('sbs1 keeps any C>T at an NCG context; sbs16 keeps any T>C at an ATN context\n')
    stop('usage: subset_by_mutsig.R in.rda {sbs1|sbs16} out.rda')
}


in.rda <- args[1]
sigtype <- args[2]
out.rda <- args[3]

if (sigtype != 'sbs1' & sigtype != 'sbs16')
    stop("sigtype must be either 'sbs1' or 'sbs16', case sensitive")

for (f in c(out.rda)) {
    if (file.exists(f))
        stop(paste('output file', f, 'already exists, please delete it first'))
}

suppressMessages(library(scan2))

if (sigtype == 'sbs1') {
    grep.pattern <- 'CG:C>T'
} else if (sigtype == 'sbs16') { 
    grep.pattern <- 'AT.:T>C'
} else {
    stop("sigtype must be either 'sbs1' or 'sbs16', case sensitive")
}

o.name <- load(in.rda)
original <- get(o.name)


if ('data.table' %in% class(original) | 'data.frame' %in% class(original)) {
    cat('detected mutation table\n')
    subsetted <- original[grep(grep.pattern, original$mutsig),]
} else if ('GenomicRanges' %in% class(original) | 'CompressedGRangesList' %in% class(original)) {
    cat('detected permutation list\n')
    subsetted <- GenomicRanges::GRangesList(lapply(original, function(x) x[grep(grep.pattern, x$mutsig),]))
} else {
    stop(paste0('object "', o.name, '" of class ', class(original), ' is not one of the recognized class types (data.table, data.frame, GenomicRanges, CompressedGRangesList)'))
}

# Saves the object named o.name, not the object o.name itself
assign(o.name, subsetted)
save(list=o.name, file=out.rda)

if ('snakemake' %in% ls()) {
    sink()
}
