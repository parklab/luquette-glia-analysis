#!/usr/bin/env Rscript

# detect script being run by snakemake
# if so, make a mock commandArgs function
if ('snakemake' %in% ls()) {
    logfile <- snakemake@log[[1]]
    con <- file(logfile, 'w')
    sink(con, type='output')
    sink(con, type='message')

    commandArgs <- function(...) unlist(c(
        snakemake@output['csv'],
        snakemake@input['muts'],
        snakemake@input['tiles'],
        snakemake@input['qbeds']
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 4) {
    cat('PASS A mutations should be used here because permutations are NOT used to control for signature bias of rescue\n')
    cat("An improved version of this tool would accept a permutation.rda file and compare correlation of actual muts to permutations, but then the resulting value wouldn't be as easy to interpret as correlation.\n")
    stop('usage: qbed_correlation.R output.csv muts.csv genome_tiles.bed qbed1 [ qbed2 ... qbedN ]')
}


out.csv <- args[1]
muts <- args[2]
tile.file <- args[3]
signal.fs <- args[-(1:3)]

for (f in c(out.csv)) {
    if (file.exists(f))
        stop(paste('output file', f, 'already exists, please delete it first'))
}

suppressMessages(library(GenomicRanges))
suppressMessages(library(data.table))
suppressMessages(library(mutenrich))

muts <- fread(muts)

count.muts <- function(tiles, x)
    countOverlaps(tiles, gr(x, add.chr.prefix=TRUE))

gr2 <- function (bed, seqinfo = NULL, add.chr.prefix = FALSE) {
    ret <- GenomicRanges::GRanges(seqnames = bed[[1]],
        ranges = IRanges::IRanges(start = bed[[2]], bed[[3]]))
    ret$keep <- bed[,5] != 0
    ret$mean.dp <- bed[[6]]
    ret
}

tiles <- gr2(fread(tile.file))

ct <- count.muts(tiles, muts)

signal.mat <- sapply(signal.fs, function(f) fread(f, skip=1)[[5]])
colnames(signal.mat) <- unname(sapply(signal.fs, function(f) mutenrich::read.bed.metadata(f, return.line=TRUE, is.qbed=TRUE)))


# Get correlation, R^2, p-values for cor=0 t-tests
values <- do.call(rbind, lapply(1:ncol(signal.mat), function(colidx) {
    col <- signal.mat[,colidx]
    df <- data.frame(signal=col[tiles$keep], mutations=ct[tiles$keep])
    m <- summary(lm(signal ~ mutations, data=df))
    data.frame(Signal=colnames(signal.mat)[colidx],
        Correlation=cor(df$signal, df$mutations),
        R.squared=m$r.squared,
        P.value=coef(m)['mutations',4])
}))

fwrite(values, file=out.csv)

if ('snakemake' %in% ls()) {
    sink()
}
