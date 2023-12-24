#!/usr/bin/env Rscript

# detect script being run by snakemake
# if so, make a mock commandArgs function
if ('snakemake' %in% ls()) {
    logfile <- snakemake@log[[1]]
    con <- file(logfile, 'w')
    sink(con, type='output')
    sink(con, type='message')

    commandArgs <- function(...) unlist(c(
        snakemake@output['pdf'],
        snakemake@output['csv'],
        snakemake@params['donor_id'],
        snakemake@input['unfilt'],
        snakemake@input[['muts']]
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}


args <- commandArgs(trailingOnly=T)

if (length(args) < 5) {
    stop("usage: ./respfig3_number_shared_snvs.R out.pdf out.csv donor_ID all___UNFILTERED_mut___any.csv shared_muts1.csv [ shared_muts2.csv ... shared_mutsN.csv ]")
}

out.pdf <- args[1]
out.csv <- args[2]
donor.id <- args[3]
unfilt.file <- args[4]
mut.files <- args[-(1:4)]

for (f in c(out.pdf, out.csv)) {
    if (file.exists(f))
        stop(sprintf("output file %s already exists, please delete it first", f))
}

# DO WORK
library(data.table)

# Not necessary to filter down to this donor
unfilt <- fread(unfilt.file) #[donor == donor.id]

muts <- lapply(mut.files, function(f) {
    print(f)
    x <- fread(f)
    x[, id := paste(chr, pos, refnt, altnt)]
    x <- unfilt[x,,on=.(id)]  # join to get filter status
    # each shared mut should have >= 2 rows, just get one representative for counting
    x[shared == TRUE, .SD[1,], by=id]
})

# Count muts passing recurrence filters
ns <- sapply(muts, function(mut) sum(!mut$final.filter))

# Get sample names
# Syntax: shared_mutations/sigprofilerextractor//5657-Oligo-3___vs___5657-Oligo-5___mutations.csv
#                                    sample 1 ---^^^^^^^^^^^^        ^^^^^^^^^^^^--- sample 2
base.names <- sapply(strsplit(mut.files, split='/'), tail, 1)  # remove preceeding directory names
sample1 <- sapply(strsplit(base.names, split='___'), `[`, 1)
sample2 <- sapply(strsplit(base.names, split='___'), `[`, 3)
fwrite(data.table(sample1=sample1, sample2=sample2, n.shared=ns), file=out.csv)

# give 1/4 em of width per pair to plot  plus a base size of 1 inch
pdf(file=out.pdf, width=1+length(mut.files)*strwidth('M', unit='inches')/4, height=2, pointsize=8)
par(mar=c(2,4,2,1))
plot(sort(ns), type='h', bty='l', xlab='',
    ylab='Number of shared sSNVs', main=paste("Subject", donor.id), xaxs='i')
points(sort(ns), pch=20, cex=3/4)
dev.off()


if ('snakemake' %in% ls()) {
    sink()
}
