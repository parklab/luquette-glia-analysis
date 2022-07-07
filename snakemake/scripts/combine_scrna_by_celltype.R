#!/usr/bin/env Rscript
#SBATCH -p park
#SBATCH -A park_contrib
#SBATCH -t 1:00:00
#SBATCH --mem=16G

# detect script being run by snakemake
# if so, make a mock commandArgs function
if ('snakemake' %in% ls()) {
    logfile <- snakemake@log[[1]]
    con <- file(logfile, 'w')
    sink(con, type='output')
    sink(con, type='message')

    commandArgs <- function(...) unlist(c(
        snakemake@input['gct'],
        snakemake@output['gct']
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}

args <- commandArgs(trailingOnly=TRUE)

if (length(args) != 2) {
    stop('usage: combine_scrnaseq_by_celltype.R scrnaseq.gct out.gct')
}

scrnaseq.gct <- args[1]
out.gct <- args[2]

if (file.exists(out.gct))
    stop(paste('output file', out.gct, 'already exists, please delete it first'))

suppressMessages(library(data.table))

rna <- as.data.frame(fread(scrnaseq.gct))

# Column names are formatted SUBJECTID.SORTTYPE___CELLTYPE
# Columns 1-2 are Ensembl gene ID, gene symbol
celltypes <- sapply(strsplit(colnames(rna)[-(1:2)], '___'), tail, 1)

out.matrix <- data.frame(rna[,c('Name', 'Description')],
    sapply(unique(celltypes), function(ct) {
        m <- as.matrix(rna[,2+which(celltypes == ct)])
        # Before summing, normalize so that libraries with more reads don't
        # overpower other libraries.
        m <- t(t(m) / colSums(m))
        rowSums(m)
    })
)
print(str(out.matrix))
colnames(out.matrix)[-(1:2)] <- unique(celltypes)

# mimic the .gct file's 2 line header
header <- c('#1.2', paste(nrow(out.matrix), ncol(out.matrix)-2, sep='\t'))

writeLines(header, con=out.gct)
write.table(out.matrix, file=out.gct,
    append=TRUE, sep='\t', col.names=TRUE, quote=F, row.names=F)

if ('snakemake' %in% ls()) {
    sink()
}
