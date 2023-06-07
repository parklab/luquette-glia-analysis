#!/usr/bin/env Rscript

# detect script being run by snakemake
# if so, make a mock commandArgs function
if ('snakemake' %in% ls()) {
    logfile <- snakemake@log[[1]]
    con <- file(logfile, 'w')
    sink(con, type='output')
    sink(con, type='message')

    commandArgs <- function(...) unlist(c(
        snakemake@input['rda'], snakemake@params['qualtype'], snakemake@output['vcf']
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 3) {
    cat('assumes in.rda contains only one object\n')
    stop("usage: scan2_to_vcf.R scan2_object.rda {A|AB|indel_A|indel_AB} out.vcf")
}

inrda <- args[1]
qualtype <- args[2]
outvcf <- args[3]

if (!(qualtype %in% c('A', 'AB', 'indel_A', 'indel_AB')))
    stop('qualtype must be one of A, AB, indel_A or indel_AB (case sensitive)')

if (file.exists(outvcf))
    stop(paste('output file', outvcf, 'already exists, please delete it first'))

suppressMessages(library(scan2))

object <- get(load(inrda))
print(object@single.cell)

# outdated name scheme for mutation calls: A = VAF-based (pass), B = mutsig rescued (rescue)
# indel_ for indels; no prefix means SNV
mt <- ifelse(substr(qualtype, 1, 5) == 'indel', 'indel', 'snv')
print( substr(qualtype, nchar(qualtype), nchar(qualtype)) )
x <- object@gatk[muttype == mt & pass == TRUE]

add.rescue <- substr(qualtype, nchar(qualtype), nchar(qualtype)) == 'B'
if (add.rescue) {
    if (is.null(object@mutsig.rescue))
        stop("requested mutsig rescued mutations, but this SCAN2 object has not undergone mutation signature based rescue")
    else
        cat("adding rescued mutations..\n")
    x <- rbind(x, object@gatk[muttype == mt & rescue == TRUE])
}


f <- file(outvcf)
# dummy header
lines.to.write <- c("##fileformat=VCFv4.0", "##source=rda_to_vcf.R", 
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">")
lines.to.write <- c(lines.to.write, paste(c("#CHROM", "POS", "ID", 
    "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", object@single.cell), 
    collapse = "\t"))

if (nrow(x) > 0)
    lines.to.write <- c(lines.to.write, paste(x$chr, x$pos, ".", x$refnt, x$altnt, ".", "PASS", ".", 'GT', '0/1', sep='\t'))

writeLines(lines.to.write, con=outvcf)

if ('snakemake' %in% ls()) {
    sink()
}
