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
        snakemake@input['scrnaseq'],
        snakemake@input['gct'],
        snakemake@output['gct']
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}

args <- commandArgs(trailingOnly=TRUE)

if (length(args) != 3) {
    stop('usage: map_scrnaseq_to_gtex_gene_model.R scrnaseq.csv gtex_tissue_median_tpm.gct out.gct')
}

scrnaseq.csv <- args[1]
gtex.median.tpm.file <- args[2]
out.gct <- args[3]

if (file.exists(out.gct))
    stop(paste('output file', out.gct, 'already exists, please delete it first'))

suppressMessages(library(data.table))

# original copy to keep order (setkey destroys it)
g <- fread(gtex.median.tpm.file, skip=2)[,1:2]
s <- fread(scrnaseq.csv)
# Join puts column 1 of g on the far right side; we want it as the first column
new.mat <- s[g, on=c('gene'='Description')]
setcolorder(new.mat, c(ncol(new.mat),1:(ncol(new.mat)-1)))
new.mat[is.na(new.mat)] <- 0

# mimic the .gct file's 2 line header
header <- c('#1.2', paste(nrow(new.mat), ncol(new.mat)-2, sep='\t'))

writeLines(header, con=out.gct)
write.table(new.mat, file=out.gct, append=TRUE, sep='\t', col.names=TRUE, quote=F, row.names=F)

if ('snakemake' %in% ls()) {
    sink()
}
