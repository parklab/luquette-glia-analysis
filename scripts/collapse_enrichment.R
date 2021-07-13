#!/usr/bin/env Rscript

# detect script being run by snakemake
# if so, make a mock commandArgs function
if ('snakemake' %in% ls()) {
    logfile <- snakemake@log[[1]]
    con <- file(logfile, 'w')
    sink(con, type='output')
    sink(con, type='message')

    commandArgs <- function(...) unlist(c(
        snakemake@input[1], snakemake@output[1:2], snakemake@params
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 4)
    stop("usage: collapse_enrichment in.rda full_output.rda summary_output.rda newclassA=oldclassA1,oldclassA2 [newclassB=oldclassB1,oldclassB2 ... ]")

inrda<- args[1]
fulloutrda <- args[2]
summaryoutrda <- args[3]
classes.to.collapse.args <- args[-(1:3)] 

classes.to.collapse <- do.call(rbind, lapply(classes.to.collapse.args, function(s) {
    x <- strsplit(s, '=')[[1]]
    newclass <- x[1]
    oldclasses <- strsplit(x[2], ',')[[1]]
    cbind(newclass, oldclasses)
}))

classes.to.collapse <- split(classes.to.collapse[,2], classes.to.collapse[,1])
print(classes.to.collapse)


if (file.exists(fulloutrda))
    stop(paste('outut file', fulloutrda, 'already exists, please delete it first'))
if (file.exists(summaryoutrda))
    stop(paste('outut file', summaryoutrda, 'already exists, please delete it first'))

source("/n/data1/hms/dbmi/park/jluquette/glia/scripts/enrichment.r")
suppressMessages(library(GenomicRanges))

load(inrda) # loads 'e'

cat('Collapsing..\n')
e <- collapse(e, classes.to.collapse=classes.to.collapse)
cat('Computing summary..\n')
es <- esummary(e)

cat('Saving..\n')
save(e, es, file=fulloutrda)
save(es, file=summaryoutrda)

if ('snakemake' %in% ls()) {
    sink()
}
