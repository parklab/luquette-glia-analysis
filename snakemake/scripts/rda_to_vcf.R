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
    cat('assumes in.rda contains only one object\n')
    stop("usage: rda_to_csv.r in.rda out.vcf")
}

inrda <- args[1]
outvcf <- args[2]

if (file.exists(outvcf))
    stop(paste('output file', outvcf, 'already exists, please delete it first'))

x <- get(load(inrda))

f <- file(outvcf)

# dummy header
vcf.header <- c("##fileformat=VCFv4.0", "##source=rda_to_vcf.R", 
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">")
vcf.header <- c(vcf.header, paste(c("#CHROM", "POS", "ID", 
    "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", 'dummy'), 
    collapse = "\t"))
writeLines(vcf.header, con = f)

writeLines(paste(x$chr, x$pos, ".", x$refnt, x$altnt, ".", "PASS", ".", 'GT', '0/1',
    sep='\t'), con=f)

if ('snakemake' %in% ls()) {
    sink()
}
