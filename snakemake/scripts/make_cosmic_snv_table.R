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
    stop("usage: make_cosmic_snv_table.R in_cosmic.csv out.csv")
}


cosmic.csv <- args[1]
outcsv <- args[2]

if (file.exists(outcsv))
    stop(paste('output file', outcsv, 'already exists, please delete it first'))

library(scan2) # for data.table and the PTA artifact signature

cosmic <- read.csv(cosmic.csv, header=T, stringsAsFactors=F)
rownames(cosmic) <- paste(cosmic[,2], cosmic[,1], sep=':')

# These are the columns that have " -    " values in them.
str(cosmic)
print(cosmic$SBS25)
print(cosmic$SBS27)
print(cosmic$SBS46)
print(cosmic$SBS47)
print(cosmic$SBS84)
print(cosmic$SBS88)
# A few columns have " -   " values for some trinuc contexts. Not sure why.
# By looking at these signatures on the COSMIC website, it seems these
# values just mean 0.
cosmic[cosmic == " -   "] <- 0
tmp.rownames <- rownames(cosmic)
cosmic <- as.data.frame(lapply(cosmic[-(1:2)], function(column) as.numeric(column)),
    stringsAsFactors=FALSE)
rownames(cosmic) <- tmp.rownames
head(cosmic)

# reordering for signatures
bases <- c("A", "C", "G", "T")
o <- paste0(rep(bases, each = 4), rep(c("C", "T"), 
        each = 48), rep(bases, times = 4), ":", rep(c("C", "T"), 
        each = 48), ">", c(rep(c("A", "G", "T"), each = 16), 
        rep(c("A", "C", "G"), each = 16)))
cosmic <- cosmic[o,]

# Add the PTA error signature
lysis.sig <- t(t(get(data(snv.artifact.signature.v3))))

if (!all(lysis.sig[,1] == rownames(cosmic)))
    stop("PTA artifact signature and COSMIC mutation types are not in the same order")

cosmic <- cbind(MutType=rownames(cosmic), cosmic, lysis.sig)

fwrite(cosmic, file=outcsv)

if ('snakemake' %in% ls()) {
    sink()
}
