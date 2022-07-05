#!/usr/bin/env Rscript
#SBATCH -p park
#SBATCH -A park_contrib
#SBATCH -t 1:00:00
#SBATCH --mem=8G

# detect script being run by snakemake
# if so, make a mock commandArgs function
if ('snakemake' %in% ls()) {
    logfile <- snakemake@log[[1]]
    con <- file(logfile, 'w')
    sink(con, type='output')
    sink(con, type='message')

    commandArgs <- function(...) unlist(c(
        snakemake@input['csv'],
        snakemake@params['sample'],
        snakemake@output['bigwig']
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}

args <- commandArgs(trailingOnly=TRUE)

if (length(args) != 3) {
    stop('usage: make_hsnp_sensitivity_track.r training_hsnps.rda sample_name out.bigwig')
}

in.rda <- args[1]
sample.name <- args[2]
out.bigwig <- args[3]

tile.size=500

if (file.exists(out.bigwig))
    stop(paste('output file', out.bigwig, 'already exists, please delete it first'))


suppressMessages(library(GenomicRanges))
suppressMessages(library(rtracklayer))
suppressMessages(library(data.table))
suppressMessages(library(mutenrich))
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg19))


cat('Reading hsnp file:', in.rda, '\n')
h <- get(load(in.rda))

cat('Got', nrow(h), 'hSNPs from sample', sample.name,
    'to build sensitivity signal for sample.\n')

cat('Tiling hg19 with', tile.size,  'bp bins\n')
bins <- tileGenome(seqlengths=seqlengths(BSgenome.Hsapiens.UCSC.hg19),
    tilewidth=tile.size, cut.last.tile.in.chrom=T)

gh <- gr(h, seqinfo=seqinfo(bins), add.chr.prefix=TRUE)
gh$pass <- TRUE  # h$pass

cat(paste0('Mapping hSNPs to non-overlapping ', tile.size, ' bp bins\n'))
snps <- countOverlaps(bins, gh)

# Make a fresh tiling GRanges with no metadata and attach the expression values for each tile.
new.bins <- GRanges(seqnames=seqnames(bins), ranges=ranges(bins), seqinfo=seqinfo(bins))
new.bins$score <- snps

cat("Writing bigwig ", out.bigwig, '\n')
export.bw(new.bins, con=out.bigwig)

if ('snakemake' %in% ls()) {
    sink()
}
