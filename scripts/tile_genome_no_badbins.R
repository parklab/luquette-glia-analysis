#!/usr/bin/env Rscript

# detect script being run by snakemake
# if so, make a mock commandArgs function
if ('snakemake' %in% ls()) {
    logfile <- snakemake@log[[1]]
    con <- file(logfile, 'w')
    sink(con, type='output')
    sink(con, type='message')

    commandArgs <- function(...) unlist(c(
        snakemake@params[1],    # binsize
        snakemake@output[1],    # out.bed
        snakemake@input['metadata'],
        snakemake@input['rdas'] # variable length list of per-chromosome RDA digest files
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 4) {
    stop("usage: tile_genome_no_badbins.R binsize out.bed metadata.yaml chr_digest1.rda [ chr.digest2.rda .. chr.digestN.rda]")
}

binsize <- as.integer(args[1])
out.bed <- args[2]
metayaml <- args[3]
digest.files <- args[-(1:3)]


if (file.exists(out.bed))
    stop(paste('output file', out.bed, 'already exists, please delete it first'))

suppressMessages(library(yaml))

meta <- read_yaml(metayaml)
sample.ids <- c(names(meta[['neuron']]), names(meta[['oligo']]))

suppressMessages(library(GenomicRanges))
suppressMessages(library(BSgenome))
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg19))

base.tiles <- do.call(c, lapply(digest.files, function(df) {
    cat('reading digest file', df, '\n')
    load(df)
    mean.mat <- mean.mat[, colnames(mean.mat) %in% sample.ids]
    tiles$dp <- rowMeans(mean.mat)
    tiles.not.in.gatk$dp <- NA
    a <- c(tiles, tiles.not.in.gatk)
    a <- sort(a)
    seqlevels(a) <- paste0('chr', seqlevels(a))
    genome(a) <- genome(BSgenome.Hsapiens.UCSC.hg19)
    a$dpclass <- ifelse(is.na(a$dp), 0,
        ifelse(a$dp < 6, 1,
            ifelse(a$dp > quantile(a$dp, prob=0.975, na.rm=T), 2, 3)))
    a
}))
seqlevels(base.tiles) <- sortSeqlevels(seqlevels(base.tiles))
#seqinfo(base.tiles) <- seqinfo(BSgenome.Hsapiens.UCSC.hg19)
seqinfo(base.tiles)

new.tiles <- tileGenome(seqlengths=seqlengths(base.tiles), #[1:22],
    tilewidth=binsize, cut.last.tile.in.chrom=T)

# Score the amount of overlap with "good" bins (class=3)
# minoverlap=100 ensures the entire bin overlaps, so that
# the number of overlapping bins can easily be interpreted as
# %overlap >= n*binsize
ols <- countOverlaps(new.tiles,
    base.tiles[base.tiles$dpclass==3,],
    minoverlap=100)   # base tiles are 100bp in size

# Require >80% overlap with passing base tiles to be included
new.tiles$status <- ifelse(ols >= 8, 'included', 'excluded')

write.table(cbind(as.character(seqnames(new.tiles)),
        start(new.tiles),
        end(new.tiles),
        1:length(new.tiles),
        ifelse(new.tiles$status=='included', 1, 0)),
    file=out.bed, quote=F, row.names=F, col.names=F, sep='\t')

if ('snakemake' %in% ls()) {
    sink()
}
