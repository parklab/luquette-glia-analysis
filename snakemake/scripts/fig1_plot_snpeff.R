#!/usr/bin/env Rscript

# detect script being run by snakemake
# if so, make a mock commandArgs function
if ('snakemake' %in% ls()) {
    logfile <- snakemake@log[[1]]
    con <- file(logfile, 'w')
    sink(con, type='output')
    sink(con, type='message')

    commandArgs <- function(...) unlist(c(
        snakemake@input[1:4], snakemake@output[1:3]
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 7) {
    stop("usage: plot_snpeff.R neuron_snpeff_SNV.vcf neuron_snpeff_INDEL.vcf oligo_snpeff_SNV.vcf oligo_snpeff_INDEL.vcf out.svg out.pdf out.csv")
}

nfile <- args[1]
nifile <- args[2]
gfile <- args[3]
gifile <- args[4]
outsvg <- args[5]
outpdf <- args[6]
outcsv <- args[7]

if (file.exists(outsvg))
    stop(paste('output file', outsvg, 'already exists, please delete it first'))
if (file.exists(outpdf))
    stop(paste('output file', outpdf, 'already exists, please delete it first'))
if (file.exists(outcsv))
    stop(paste('output file', outcsv, 'already exists, please delete it first'))

suppressMessages(library(data.table))
suppressMessages(library(svglite))

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

res <- rbind(
    get.severity(nfile),
    get.severity(nifile),
    get.severity(gfile),
    get.severity(gifile))
df <- data.frame(CellType=c("Neuron", "Neuron", 'Oligo', 'Oligo'),
    MutType=c('SNV', 'Indel', 'SNV', 'Indel'),
    res, stringsAsFactors=F)
totals <- rowSums(res)
    

figwidth=5
figheight=3
devs <- list(svglite, pdf)
outs <- c(outsvg, outpdf)
for (i in 1:2) {
    devs[[i]](file=outs[i], width=figwidth, height=figheight)

    layout(t(1:2))
    par(mar=c(6,4,3,1))
    snvs <- 100*res[c(1,3),-4]/totals[c(1,3)]
    indels <- 100*res[1+c(1,3),-4]/totals[1+c(1,3)]
    ylim=c(0, max(c(snvs, indels)))
    barplot(snvs, ylim=ylim, beside=T,
        border=NA, col=1:2, main='SNV', ylab='Percent', las=3)
    barplot(indels, ylim=ylim, beside=T,
        border=NA, col=1:2, main='Indel', ylab='Percent', las=3)
}

fwrite(df, file=outcsv)

if ('snakemake' %in% ls()) {
    sink()
}
