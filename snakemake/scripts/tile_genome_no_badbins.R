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
    stop("usage: tile_genome_no_badbins.R binsize out.bed metadata.csv chr_digest1.rda [ chr.digest2.rda .. chr.digestN.rda]")
}

binsize <- as.integer(args[1])
out.bed <- args[2]
metacsv <- args[3]
digest.files <- args[-(1:3)]


if (file.exists(out.bed))
    stop(paste('output file', out.bed, 'already exists, please delete it first'))

suppressMessages(library(data.table))

meta <- fread(metacsv)
sample.ids <- meta$sample
print(sample.ids)

suppressMessages(library(GenomicRanges))
suppressMessages(library(BSgenome))
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg19))

# unname: IMPORTANT - snakemake names the digest files, which, for reasons
# I don't understand, prevents GRanges from c()ing together.
base.tiles <- do.call(c, lapply(unname(digest.files), function(df) {
    cat('reading digest file', df, '\n')
    load(df)
    mean.mat <- mean.mat[, ..sample.ids]
    tiles$dp <- rowMeans(mean.mat)
    tiles.not.in.gatk$dp <- NA
    a <- c(tiles, tiles.not.in.gatk)
    a <- sort(a)
    seqlevels(a) <- paste0('chr', seqlevels(a))
    seqlevels(a) <- seqlevels(BSgenome.Hsapiens.UCSC.hg19)
    seqinfo(a) <- seqinfo(BSgenome.Hsapiens.UCSC.hg19)
    a$dpclass <- ifelse(is.na(a$dp), 0,
        ifelse(a$dp < 6, 1,
            ifelse(a$dp > quantile(a$dp, prob=0.975, na.rm=T), 2, 3)))
    a$dp[is.na(a$dp)] <- 0   # for mean computations later, better to have 0. the NA status is already saved in class=0
    a
}))

new.tiles <- tileGenome(seqlengths=seqlengths(base.tiles)[paste0('chr', 1:22)],
    tilewidth=binsize, cut.last.tile.in.chrom=T)

# Score the amount of overlap with "good" bins (class=3)
# minoverlap=100 ensures the entire bin overlaps, so that
# the number of overlapping bins can easily be interpreted as
# %overlap >= n*binsize
cat('counting overlaps..\n')
nols <- countOverlaps(new.tiles,
    base.tiles[base.tiles$dpclass==3,],
    minoverlap=100)   # base tiles are 100bp in size

# Require >80% overlap with passing base tiles to be included
new.tiles$status <- ifelse(nols >= 0.8*binsize/100, 'included', 'excluded')


# Get mean depth for each tile. In this case we don't consider only dpclass==3
# base tiles.
cat('finding overlaps for mean depth..\n')
ols <- findOverlaps(new.tiles, base.tiles, minoverlap=100)
ols <- split(to(ols), from(ols))
new.tiles$mean.dp <- 0   # anything not in base.tiles treated as 0 depth
new.tiles[as.integer(names(ols)),]$mean.dp <- lapply(ols, function(tos) mean(base.tiles[tos,]$dp))

write.table(cbind(as.character(seqnames(new.tiles)),
        start(new.tiles),
        end(new.tiles),
        1:length(new.tiles),
        ifelse(new.tiles$status=='included', 1, 0),
        new.tiles$mean.dp),
    file=out.bed, quote=F, row.names=F, col.names=F, sep='\t')

if ('snakemake' %in% ls()) {
    sink()
}
