#!/usr/bin/env Rscript

# detect script being run by snakemake
# if so, make a mock commandArgs function
if ('snakemake' %in% ls()) {
    logfile <- snakemake@log[[1]]
    con <- file(logfile, 'w')
    sink(con, type='output')
    sink(con, type='message')

    commandArgs <- function(...) unlist(c(
        snakemake@input['rds'],
        snakemake@output['csv']
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}


args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 2) {
    cat("This script requires the seurat environment.\n")
    cat("Correction: if we manually access the Seurat object (avoid using Seurat methods), only the data.table library is necessary for this script.\n")
    stop('usage: rds_to_csv.R in.rds out.csv')
}


rds.file <- args[1]
out.csv <- args[2]

if (file.exists(out.csv))
    stop(paste('output file', out.csv, 'already exists, please delete it first'))

#library(Seurat)
library(data.table)

# x is a Seurat object
x <- readRDS(rds.file)

# Access slots directly to avoid loading library(Seurat), do not use any methods (e.g., Idents())
z <- data.table(Library=x@meta.data$samples,
    #CellType=as.character(Idents(x)),
    CellType=as.character(x@meta.data$CellTypes),
    x@reductions$umap@cell.embeddings)

fwrite(z, file=out.csv)

if ('snakemake' %in% ls()) {
    sink()
}
