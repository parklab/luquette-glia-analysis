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
if (length(args) < 6) {
    cat('hg19 is assumed; chromosome should be prefixed with "chr"\n')
    cat('base.tile.size is the smallest tile size used to create a\n')
    cat('non-overlapping tiling of the genome. Windows created further\n')
    cat('down in the pipeline CAN ONLY be integer multiples of this tile\n')
    cat('size. So the tile size here should be very small; e.g., 100 bp.\n')
    stop("usage: digest_chrom_depth.R chromosome base.tile.size tiles.per.chunk n.cores out.rda joint_depth_matrix1.tab.gz [ joint_depth_matrix2.tab.gz ... joint_depth_matrixN.tab.gz ]")
}

chrom <- args[1]
base.tile.size <- as.integer(args[2])
tiles.per.chunk <- as.integer(args[3])
n.cores <- as.integer(args[4])
out.rda <- args[5]
matfiles <- args[-(1:5)]

if (file.exists(out.rda))
    stop(paste('output file', out.rda, 'already exists, please delete it first'))

suppressMessages(library(scan2))
suppressMessages(library(pryr))
suppressMessages(library(progressr))
suppressMessages(library(future))
suppressMessages(library(future.apply))
suppressMessages(library(GenomeInfoDb))
if (n.cores > 1)
    plan(multicore, workers=n.cores)

genome <- GenomeInfoDb::Seqinfo(genome='hg19')

if (!(chrom %in% seqlevels(genome)))
    stop(paste('chromosome', chrom, 'is not a valid hg19 chromosome. Did you add the "chr" prefix?'))

# Promote warnings to errors in (for catching GenomicRanges warnings)
options(warn=2)

# First step: build the full tile map over the chromosome
# Build a base set of tiles covering the genome at the maximum
# resolution (base.tile.size).
# autosomes only
chrom.end <- seqlengths(genome)[chrom]
chrom.end <- min(chrom.end, 1e6)

tiles <- tileGenome(seqlengths=setNames(chrom.end, chrom),
    tilewidth=base.tile.size, cut.last.tile.in.chrom=T)


# build a GRanges object of positions from GATK DepthOfCoverage.
# GATK DepthOfCoverage outputs many 0 values, but it skips over
# some parts of the genome (probably Ns).
positions <- read.tabix.data(path=matfiles[1],
    region=GRanges(seqnames=sub('chr', '', chrom),  # SCAN2 output doesn't have chr prefix
                   ranges=IRanges(start=1, end=chrom.end)),
    colClasses=c('NULL', 'integer'))[[1]]
g.basepair <- GRanges(seqnames=chrom, ranges=IRanges(start=positions, width=1))

# There are a handful of cases where 1-5 bases are not reported
# by GATK. These shouldn't cause any problems, so use min.gapwidth to
# merge over them.
gatk <- reduce(g.basepair, min.gapwidth=50)


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
# the 95% overlap necessary to keep the tile. This affects very little
# of the genome due to base.tile.size being small.
tiles.not.in.gatk <- subsetByOverlaps(tiles, gatk, minoverlap=0.95*base.tile.size,
    invert=TRUE)
tiles <- subsetByOverlaps(tiles, gatk, minoverlap=0.95*base.tile.size)
print(tiles)


# Create a mapping from g.basepair -> tiles. This only needs to be
# done once since all matfiles have the same position format.
cat('Creating position -> tile map\n')
system.time(tilemap <- findOverlaps(g.basepair, tiles))
str(tilemap)
cat('Tilemap size in MB (must be duplicated for each thread):', object_size(tilemap)/1e6, '\n')


cat('Memory profile before matrix construction:\n')
print(gc())

tiles$chunk.id <- head(rep(1:ceiling(length(tiles)/tiles.per.chunk), each=tiles.per.chunk), length(tiles))
chunks <- unlist(reduce(split(tiles, tiles$chunk.id)))
# SCAN2 files don't have chr prefixes
seqnames(chunks) <- sub('chr', '', seqnames(chunks))
cat('Starting depth matrix digestion on', length(chunks), 'chunks.\n')
cat('Parallelizing using', future::nbrOfWorkers(), 'cores.\n')
cat('chunks:\n')
print(chunks)
cat("files:\n")
print(matfiles)

# IMPORTANT!
# tilemap has to be exported to child processes. This can be very large
# if tiles.per.chunk is too big.
progressr::with_progress({
    p <- progressr::progressor(along=1:(length(chunks)*length(matfiles)))
    p(amount=0, class='sticky', scan2::perfcheck(print.header=TRUE))
    mean.mat <- rbindlist(future.apply::future_lapply(1:length(chunks), function(i) {
        # Create a matrix of average read depth in each tile for each sample.
        mean.mats <- lapply(1:length(matfiles), function(j) {
            pc <- perfcheck(paste('chunk', i, 'file', j, '/', length(matfiles)), {
                f <- matfiles[j]
                dpm.basepair <- read.tabix.data(f, region=chunks[i])
                if (nrow(dpm.basepair) == 0) {
                    ret <- c()
                } else {
                    # initialize to NA because all bp positions are not in the tile map
                    dpm.basepair[, tileid := NA]  
                    dpm.basepair[from(tilemap), tileid := to(tilemap)]
                    ret <- dpm.basepair[,setNames(as.list(colMeans(.SD)), colnames(.SD)),by=tileid][!is.na(tileid)][, -'tileid']
                }
            })
            p(class='sticky', amount=1, pc)
            ret
        })
        do.call(cbind, mean.mats)
    }))
}, enable=TRUE)

save(chrom, base.tile.size, tiles, tiles.not.in.gatk, mean.mat,
    file=out.rda, compress=FALSE)

cat('Final memory profile:\n')
print(gc())

if ('snakemake' %in% ls()) {
    sink()
}
