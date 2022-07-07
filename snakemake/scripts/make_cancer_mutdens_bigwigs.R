#!/usr/bin/env Rscript

#SBATCH -p short
#SBATCH -t 0:30:00
#SBATCH --mem=32G
#SBATCH -c 1

# detect script being run by snakemake
# if so, make a mock commandArgs function
if ('snakemake' %in% ls()) {
    logfile <- snakemake@log[[1]]
    con <- file(logfile, 'w')
    sink(con, type='output')
    sink(con, type='message')

    commandArgs <- function(...) unlist(c(
        snakemake@input['txt'],
        snakemake@output['sumbigwig'],
        snakemake@output['normbigwig']
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 3) {
    stop('[script.r] cancer_mut.rda out_sumdens.bigWig out_normdens.bigWig')
}


#cancer.rda <- args[1]
cancer.txt <- args[1]
out.sumdens <- args[2]
out.normdens <- args[3]

if (file.exists(out.sumdens))
    stop(paste('output file', out.sumdens, 'already exists, please delete it first'))
if (file.exists(out.normdens))
    stop(paste('output file', out.normdens, 'already exists, please delete it first'))


suppressMessages(library(data.table))
suppressMessages(library(GenomicRanges))
suppressMessages(library(mutenrich))  # for gr()
suppressMessages(library(BSgenome))
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg19))

bins <- tileGenome(seqlengths=seqlengths(BSgenome.Hsapiens.UCSC.hg19)[sortSeqlevels(seqlevels(BSgenome.Hsapiens.UCSC.hg19))[1:22]],
    tilewidth=100, cut.last.tile.in.chrom=T)

count.muts <- function(tiles, x)
    countOverlaps(tiles, gr(x, add.chr.prefix=TRUE))


# convert genome tile file into GR obj
gr2 <- function (bed, seqinfo = NULL, add.chr.prefix = FALSE) {
    ret <- GenomicRanges::GRanges(seqnames = bed[[1]], ranges = IRanges::IRanges(start = bed[[2]], bed[[3]]))
    ret$keep <- bed[,5] != 0
    ret
}


summarize.muts <- function(mut, bins, normalize.by.sample=FALSE) {
    # normalize each sample's mut counts such that the sum of signal contributed over
    # all windows is 1
    if (normalize.by.sample) {
        cat('weight: sample table\n')
        print(table(mut$sample))
        ret <- sapply(split(mut, mut$sample), count.muts, tiles=bins)
        rowSums(t(t(ret)/colSums(ret)))
    } else {
        count.muts(bins, mut)
    }
}

cat('reading', cancer.txt, '\n')
mut <- fread(cancer.txt)

# determine which tumors are hypermutated by Tukey's criterion
# Remove hypermutated samples using Tukey's IQR*1.5 criteria.
cat("Outlier detection (Tukey's IQR*1.5)\n")
nmuts <- table(mut$sample)
print(sort(nmuts))
q3 <- quantile(nmuts, prob=3/4)
iqr <- IQR(nmuts)
cutoff <- q3 + 1.5*iqr
print(summary(as.vector(nmuts)))
cat('    q3 = ', q3, '\n')
cat('    IQR = ', iqr, '\n')
cat("    cutoff = ", cutoff, '\n')
cat('    removing', sum(nmuts > cutoff), 'samples\n')
outliers <- names(nmuts)[nmuts > cutoff]
mut <- mut[!(mut$sample %in% outliers),]


# mark recurrent sites; removing later
mut$recurrent <- duplicated(paste(mut$chr, mut$pos))
mut <- mut[!mut$recurrent,]


summut <- summarize.muts(mut, bins, normalize.by.sample=FALSE)
normmut <- summarize.muts(mut, bins, normalize.by.sample=TRUE)


cat("writing sum density bigwig..\n")
bins$score <- summut
print(bins)
export.bw(bins, out.sumdens)

cat("writing normalized density bigwig..\n")
bins$score <- normmut
print(bins)
export.bw(bins, out.normdens)
