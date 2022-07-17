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
        snakemake@input['background'],
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
#print('before normalization')
#print(colSums(m, na.rm=TRUE))
background.matrix <- t(t(m)/colSums(m, na.rm=TRUE))  # normalize read counts per cell type
#print('after normalization')
#print(colSums(background.matrix, na.rm=TRUE))
mean.background <- rowMeans(background.matrix, na.rm=TRUE)

# convert the score to a fold change
# use a small epsilon to avoid division by 0
eps <- min(background.matrix[background.matrix>0], na.rm=TRUE)/2
out[[5]] <- out[[5]]/sum(out[[5]], na.rm=TRUE) # not necessary for quantile based analysis, but nice to have for interpretation
out[[5]] <- (out[[5]]+eps) / (mean.background+eps)
# compute quantiles
out[[4]] <- findInterval(out[[5]],
    quantile(out[[5]], na.rm=T, probs=1:nq/nq), rightmost.closed=TRUE)+1

writeLines(ins[[target]]$metadata, out.file)
fwrite(out, file=out.file, append=TRUE, quote=FALSE, sep='\t', col.names=FALSE)

if ('snakemake' %in% ls()) {
    sink()
}
