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
        snakemake@input['sens'],
        snakemake@input['tiles'],
        snakemake@input['qbeds']
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 5) {
    cat('PASS A mutations should be used here because permutations are NOT used to control for signature bias of rescue\n')
    cat("An improved version of this tool would accept a permutation.rda file and compare correlation of actual muts to permutations, but then the resulting value wouldn't be as easy to interpret as correlation.\n")
    stop('usage: qbed_correlation.R output.csv muts.csv sens.summary.txt genome_tiles.bed qbed1 [ qbed2 ... qbedN ]')
}


out.csv <- args[1]
muts <- args[2]
sens.file <- args[3]
tile.file <- args[4]
signal.fs <- args[-(1:4)]

for (f in c(out.csv)) {
    if (file.exists(f))
        stop(paste('output file', f, 'already exists, please delete it first'))
}

suppressMessages(library(GenomicRanges))
suppressMessages(library(data.table))
suppressMessages(library(mutenrich))

muts <- fread(muts)
muttype <- muts$muttype[1]
if (any(muts$muttype != muttype))
    stop('this script can only be applied to one mutation type at a time')

if (any(muts$rescue))
    stop('this script can only be applied to VAF-based SCAN2 calls')

if (any(muts$pass == FALSE))
    stop("this script can only be applied to VAF-based SCAN2 calls. all should have pass=TRUE")

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

sens <- fread(sens.file)
# sensitivity tables are sorted lexicographically. reorder them to match tiles
sens <- sens[order(as.integer(chr), start)]

if (any(seqnames(tiles) != paste0('chr',sens$chr)))
    stop("column CHR does not match for tiles and sensitivity table")
if (any(start(tiles) != sens$start))
    stop("column START does not match for tiles and sensitivity table")
if (any(end(tiles) != sens$end))
    stop("column END does not match for tiles and sensitivity table")

ct <- count.muts(tiles, muts)
# Correct mutation count for sensitivity by redistributing the same total
# number of mutations over the bins at relative rates defined by
# weighted.<muttype>.sens
sens.column <- paste0('weighted.', muttype, '.sens')
total.muts <- nrow(muts)
if (total.muts != sum(ct))
    stop('number of intersections of mutations with tiles != number of mutations. is the tiling disjoint?')

print(total.muts)
s <- sens[[sens.column]]
print(s)
# scaling factor
uncorrected.ct <- ct
ct <- ifelse(s==0,0,ct/s) * sum(ifelse(s==0,0,ct))/sum(ifelse(s==0,0,ct/s))
print(ct)

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
