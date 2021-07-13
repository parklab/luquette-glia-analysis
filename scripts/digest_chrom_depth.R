#!/usr/bin/env Rscript

# detect script being run by snakemake
# if so, make a mock commandArgs function
if ('snakemake' %in% ls()) {
    logfile <- snakemake@log[[1]]
    con <- file(logfile, 'w')
    sink(con, type='output')
    sink(con, type='message')

    commandArgs <- function(...) unlist(c(
        snakemake@params[1:2],
        snakemake@output[1],
        snakemake@input[1], # position file
        snakemake@input[-1] # variable number of depth files
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 5) {
    cat('hg19 is assumed; chromosome should be prefixed with "chr"\n')
    cat('base.tile.size is the smallest tile size used to create a\n')
    cat('non-overlapping tiling of the genome. Windows created further\n')
    cat('down in the pipeline CAN ONLY be integer multiples of this tile\n')
    cat('size.\n')
    stop("usage: digest_chrom_depth.R chromosome base.tile.size out.rda positions.txt depth1.txt [ depth2.txt ... depthN.txt ]")
}

chrom <- args[1]
base.tile.size <- as.integer(args[2])
out.rda <- args[3]
posfile <- args[4]
dpfiles <- args[-(1:4)]

if (file.exists(out.rda))
    stop(paste('output file', out.rda, 'already exists, please delete it first'))

suppressMessages(library(data.table))
suppressMessages(library(GenomicRanges))
suppressMessages(library(BSgenome))
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg19))

if (!(chrom %in% seqlevels(BSgenome.Hsapiens.UCSC.hg19)))
    stop(paste('chromosome', chrom, 'is not a valid hg19 chromosome. Did you add the "chr" prefix?'))

# GenomicRanges only warns when two objects don't have compatible
# chromosome names. Make it cause an error so we fail early.
options(warn=2)

# build a GRanges object of positions from GATK DepthOfCoverage.
# GATK DepthOfCoverage outputs many 0 values, but it skips over
# some parts of the genome (probably Ns).
positions <- fread(posfile)[[1]]
g.basepair <- GRanges(seqnames=chrom, ranges=IRanges(start=positions, width=1))

# There are a handful of cases where 1-5 bases are not reported
# by GATK. These shouldn't cause any problems, so use min.gapwidth to
# merge over them.
gatk <- reduce(g.basepair, min.gapwidth=50)


# Build a base set of tiles covering the genome at the maximum
# resolution (base.tile.size).
seqinfo(BSgenome.Hsapiens.UCSC.hg19)
# autosomes only
tiles <- tileGenome(seqlengths=seqlengths(BSgenome.Hsapiens.UCSC.hg19)[chrom],
    tilewidth=base.tile.size, cut.last.tile.in.chrom=T)


# Keep any tile that overlaps at least 95% with GATK coverage.
# This doesn't quite do that, but because of the reduce() above, there
# should be very few ranges in the 'gatk' object and they should be
# quite large. This code will exclude windows that, for example, overlap
# by >95% with GATK coverage but the GATK coverage is broken up into
# two smaller windows such that each individual window doesn't give 95%
# overlap. E.g.,
#          |---------------- tile ------------------|
#   |---- gatk window 1 -------| |------ gatk window 2 -------|
# In the above example, the tile overlaps >95% with GATK, but because
# that overlap is computed per GATK window, neither windows 1 or 2 meet
# the 95% overlap necessary to keep the tile. This should happen very
# rarely.
tiles.not.in.gatk <- subsetByOverlaps(tiles, gatk, minoverlap=0.95*base.tile.size,
    invert=TRUE)
tiles <- subsetByOverlaps(tiles, gatk, minoverlap=0.95*base.tile.size)


# Create a mapping from g.basepair -> tiles. This only needs to be
# done once since all files have the same position format.
cat('Creating position -> tile map\n')
system.time(tilemap <- findOverlaps(g.basepair, tiles))


sample.ids <- sub('.txt$', '', basename(dpfiles)) # use file name as sample ID
cat('Memory profile before matrix construction:\n')
print(gc())

# Create a matrix of average read depth in each tile for each sample.
mean.dp.per.tile.matrix <- sapply(setNames(dpfiles, sample.ids), function(f) {
    cat('reading', f, '\n')
    dp.basepair <- fread(f)[[1]]  # fread is way faster than scan() for some reason
    sample.id <- sub('.txt$', '', basename(f)) # use file name as sample ID
    mean.dp.per.tile <- sapply(split(dp.basepair[from(tilemap)], to(tilemap)), mean)
    mean.dp.per.tile
})

save(chrom, base.tile.size, tiles, tiles.not.in.gatk, mean.dp.per.tile.matrix,
    file=out.rda)

cat('Final memory profile:\n')
print(gc())

if ('snakemake' %in% ls()) {
    sink()
}
