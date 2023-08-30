#!/usr/bin/env Rscript

# detect script being run by snakemake
# if so, make a mock commandArgs function
if ('snakemake' %in% ls()) {
    logfile <- snakemake@log[[1]]
    con <- file(logfile, 'w')
    sink(con, type='output')
    sink(con, type='message')

    commandArgs <- function(...) unlist(c(
        snakemake@params['group_tag'],
        snakemake@input['vcf'],
        snakemake@output['csv']
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 3) {
    cat("group.tag is added to each output row to identify which cells were used to generate snpeff.vcf\n")
    stop("usage: plot_snpeff.R group.tag snpeff.vcf out.csv")
}

group.tag <- args[1]
in.vcf <- args[2]
out.csv <- args[3]

if (file.exists(out.csv))
    stop(paste('output file', out.csv, 'already exists, please delete it first'))

suppressMessages(library(data.table))

levels <- c('HIGH', 'MODERATE', 'LOW', 'MODIFIER')

# The INFO tag in these VCFs starts with 3 tags: TType={pta_neuron|pta_oligo|...};Origin={something};Sample={something};ANN=...
#
# The element starting with ANN=... is the SnpEff string.
get.snpeff.string <- function(s)
    paste(sapply(strsplit(s, ';', fixed=TRUE), function(x) x[-(1:3)]), sep=';')
get.severity.string <- function(s)
    sapply(strsplit(s, '\\|'), function(x) x[3])
get.severity <- function(filename) {
    x <- read.table(filename, comment='#', sep='\t', stringsAsFactors=FALSE)
    y <- get.snpeff.string(x[,8])
    z <- get.severity.string(y)
    sapply(levels, function(lv) sum(z == lv, na.rm=TRUE))
}

sev <- get.severity(in.vcf)

# Mimic enrichment output so we can use barplot_enrich_new2.R
dt <- data.table(group.tag=group.tag, class=names(sev), enr=100*sev/sum(sev), enr.boot.0.95.lb=NA, enr.boot.0.95.ub=NA)

fwrite(dt[class != 'MODIFIER'], file=out.csv)

if ('snakemake' %in% ls()) {
    sink()
}
