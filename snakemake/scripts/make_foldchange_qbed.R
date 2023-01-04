#!/usr/bin/env Rscript

# detect script being run by snakemake
# if so, make a mock commandArgs function
if ('snakemake' %in% ls()) {
    logfile <- snakemake@log[[1]]
    con <- file(logfile, 'w')
    sink(con, type='output')
    sink(con, type='message')

    commandArgs <- function(...) unlist(c(
        snakemake@output['qbed'],
        snakemake@input['target'],
        snakemake@input['background']
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}

suppressMessages(library(data.table))
suppressMessages(library(mutenrich))

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 3) {
    stop("usage: make_foldchange_qbed.R output_file target_file background_file1 [ background_file2 ... background_fileN ]")
}

files <- unique(args[-1]) # just in case user supplies target as a background file
target <- args[2]
out.file <- args[1]

use.old <- FALSE

if (file.exists(out.file))
    stop(paste('output file', out.file, 'already exists, please delete it first'))

ins <- lapply(files, function(f) {
    list(metadata=readLines(f, n=1),
        data=fread(f, skip=1))
})
names(ins) <- files

out <- ins[[target]]$data
out.metadata <- read.bed.metadata(target, is.qbed=TRUE)
# add _foldchange to the datasource
out.metadata$datasource <- paste0(out.metadata$datasource, '_foldchange')
nq=as.integer(out.metadata$QUANTILES)

# don't include the target track in the background
m <- sapply(ins[names(ins) != target], function(i) i$data[[5]])

background.matrix <- t(t(m)/colSums(m, na.rm=TRUE))  # normalize read counts per cell type
mean.background <- rowMeans(background.matrix, na.rm=TRUE)

# This older method had 2 issues when dealing with scRNAseq data:
#   1. without filtering for the top/bottom 2.5% and number of zeros
#      in the input matrix, a large number of (0+eps)/(0+eps) sites
#      existed. such sites are more likely to reflect poor accessibility
#      to the assay than useful data points.
#   2. using normalized count / mean normalized count ignores the variance
#      for each bin. we now use (normalized count - mean normalized count) / sd(normalized count)
#      to better account for this. the result of ignoring variance was that
#      bins 1-2 and 9-10 had low mean normalized count compared to bins 3-8
#      (plotting them looked like an upside-down U shape). this likely meant
#      that bin assignment at the extreme ends was driven more by unaccounted-for
#      variance than true signal. after this normalization, mean normalized
#      count distribution across bins 1-10 greatly improved, though it is still
#      not as flat as one would hope for.
if (!use.old) {
    nzeros <- rowSums(background.matrix == 0)  # counts how many cell types have 0 entries
    under.2.5 <- mean.background < quantile(mean.background, na.rm=TRUE, prob=0.025)
    over.97.5 <- mean.background > quantile(mean.background, na.rm=TRUE, prob=0.975)
    bin.keep <- nzeros <= 2 & under.2.5 == FALSE & over.97.5 == FALSE
    bin.keep[is.na(bin.keep)] <- FALSE
    # get rid of bins with 3 or more 0s or with bottom 2.5%/top 2.5% of signal value.
    bm <- t(t(background.matrix)/colSums(background.matrix[bin.keep,], na.rm=TRUE))  # normalize read counts per cell type

    ms <- rowMeans(bm, na.rm=TRUE)
    sds <- apply(bm, 1, sd)
    n <- out[[5]] / sum(out[[5]][bin.keep], na.rm=TRUE)  # normalize to sum(count)=1, but only over bins passing above thresholds
    msd <- (n - ms) / sds   # quantity that accounts for both the background mean and variance across cell types
    q <- findInterval(msd, quantile(msd[bin.keep], na.rm=T, probs=1:nq/nq), rightmost.closed=T)+1
    q[!bin.keep] <- 'foldchange_removed'  # special quantile name for bins removed by above filters

    # store scores/quantiles in expected QBED columns
    out[[5]] <- msd
    out[[4]] <- q
} else if (use.old) {
    # convert the score to a fold change
    # use a small epsilon to avoid division by 0
    eps <- min(background.matrix[background.matrix>0], na.rm=TRUE)/2
    out[[5]] <- out[[5]]/sum(out[[5]], na.rm=TRUE) # not necessary for quantile based analysis, but nice to have for interpretation
    out[[5]] <- (out[[5]]+eps) / (mean.background+eps)
    # compute quantiles
    out[[4]] <- findInterval(out[[5]],
        quantile(out[[5]], na.rm=T, probs=1:nq/nq), rightmost.closed=TRUE)+1
}

writeLines(ins[[target]]$metadata, out.file)
fwrite(out, file=out.file, append=TRUE, quote=FALSE, sep='\t', col.names=FALSE)

if ('snakemake' %in% ls()) {
    sink()
}
