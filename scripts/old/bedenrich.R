#!/usr/bin/env Rscript
#SBATCH -p park
#SBATCH -A park_contrib
#SBATCH -t 1:00:00
#SBATCH --mem=10G

# detect script being run by snakemake
# if so, make a mock commandArgs function
if ('snakemake' %in% ls()) {
    logfile <- snakemake@log[[1]]
    con <- file(logfile, 'w')
    sink(con, type='output')
    sink(con, type='message')

    commandArgs <- function(...) unlist(c(
        snakemake@params['nbootstraps'],
        snakemake@input['mut'],
        snakemake@input['perm'],
        snakemake@output['full'],
        snakemake@output['summary'],
        snakemake@input['qbed']
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}


suppressMessages(library(mutenrich))

command.line.analysis(function(genome, bedf)
    enrich.data(gbed=read.bed(bedf, genome, feature.name='feature', add.chr.prefix=TRUE)),
    genome="BSgenome.Hsapiens.UCSC.hg19", #"BSgenome.Hsapiens.1000genomes.hs37d5",
    args=commandArgs(trailingOnly=TRUE))
