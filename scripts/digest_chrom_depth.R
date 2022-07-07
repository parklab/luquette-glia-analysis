#!/usr/bin/env Rscript

# detect script being run by snakemake
# if so, make a mock commandArgs function
if ('snakemake' %in% ls()) {
    logfile <- snakemake@log[[1]]
    con <- file(logfile, 'w')
    sink(con, type='output')
    sink(con, type='message')

    commandArgs <- function(...) unlist(c(
        snakemake@params['chrom'],
        snakemake@params['base_tile_size'],
        snakemake@params['chunks_per_tile'],
        snakemake@threads,
        snakemake@output['rda'],
        snakemake@input['matfiles'] # variable number of matrix files
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 6) {
    cat('hg19 is assumed.\n')
    cat('base.tile.size is the smallest tile size used to create a\n')
    cat('non-overlapping tiling of the genome. Windows created further\n')
    cat('down in the pipeline CAN ONLY be integer multiples of this tile\n')
    cat('size. So the tile size here should be very small; e.g., 100 bp.\n')
    stop("usage: digest_chrom_depth.R chromosome base.tile.size tiles.per.chunk n.cores out.rda joint_depth_matrix1.tab.gz [ joint_depth_matrix2.tab.gz ... joint_depth_matrixN.tab.gz ]")
}

chrom <- as.character(args[1])
base.tile.size <- as.integer(args[2])
tiles.per.chunk <- as.integer(args[3])
n.cores <- as.integer(args[4])
out.rda <- args[5]
matfiles <- args[-(1:5)]

if (file.exists(out.rda))
    stop(paste('output file', out.rda, 'already exists, please delete it first'))

suppressMessages(library(scan2))
suppressMessages(library(progressr))
suppressMessages(library(future))
suppressMessages(library(future.apply))
suppressMessages(library(GenomeInfoDb))
if (n.cores > 1)
    plan(multicore, workers=n.cores)

genome <- GenomeInfoDb::Seqinfo(genome='GRCh37.p13')

if (!(chrom %in% seqlevels(genome)))
    stop(paste('chromosome', chrom, 'is not a valid GRCh37.p13 chromosome'))

# Promote warnings to errors in (for catching GenomicRanges warnings)
options(warn=2)

# First step: build the full tile map over the chromosome
# Build a base set of tiles covering the genome at the maximum
# resolution (base.tile.size).
# autosomes only
chrom.end <- seqlengths(genome)[chrom]

tiles <- tileGenome(seqlengths=setNames(chrom.end, chrom),
    tilewidth=base.tile.size, cut.last.tile.in.chrom=T)


build.tilemap <- function(chunk, tiles, representative.matfile) {
    # build a GRanges object of positions from GATK DepthOfCoverage.
    # GATK DepthOfCoverage outputs many 0 values, but it skips over
    # some parts of the genome (probably Ns).
    positions <- read.tabix.data(path=representative.matfile,
        region=chunk, colClasses=c('NULL', 'integer'))#[[1]]
    if (nrow(positions) == 0)
        return(NULL)

    g.basepair <- GRanges(seqnames=seqnames(tiles)[1],
        ranges=IRanges(start=positions[[1]], width=1))
    # XXX: see discussion about tabix bug below. tabix returns positions
    # outside of the requested chunk so we have to remove them.
    g.basepair <- subsetByOverlaps(g.basepair, chunk)

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
    
    tilemap <- findOverlaps(g.basepair, tiles)
    list(tiles=tiles, tiles.not.in.gatk=tiles.not.in.gatk, tilemap=tilemap)
}


cat('Memory profile before matrix construction:\n')
print(gc())

tiles$chunk.id <- head(rep(1:ceiling(length(tiles)/tiles.per.chunk), each=tiles.per.chunk), length(tiles))
chunks <- unlist(reduce(split(tiles, tiles$chunk.id)))
cat('Starting depth matrix digestion on', length(chunks), 'chunks.\n')
cat('Parallelizing using', future::nbrOfWorkers(), 'cores.\n')
cat('chunks:\n')
print(chunks)
cat("files:\n")
print(matfiles)


progressr::with_progress({
    p <- progressr::progressor(along=1:(length(chunks)*length(matfiles)))
    progressr::handlers(progressr::handler_newline())
    p(amount=0, class='sticky', scan2::perfcheck(print.header=TRUE))
    results <- future.apply::future_lapply(1:length(chunks), function(i) {
        tm <- build.tilemap(chunk=chunks[i], tiles=subsetByOverlaps(tiles, chunks[i]),
            representative.matfile=matfiles[1])
        # returns NULL for regions with no data; e.g., heterochromatin arms like 22p
        if (is.null(tm)) {
            tiles <- c()
            tiles.not.in.gatk <- c()
            mean.mat <- c()
            cat("skipping chunk", i, "; no depth data. This is normal for heterochromatin arms like 22p.\n")
            print(chunks[i])
        } else {
            tilemap <- tm$tilemap
            tiles <- tm$tiles
            tiles.not.in.gatk <- tm$tiles.not.in.gatk
            # Create a matrix of average read depth in each tile for each sample.
            mean.mats <- lapply(1:length(matfiles), function(j) {
                pc <- perfcheck(paste('chunk', i, 'file', j, '/', length(matfiles)), {
                    f <- matfiles[j]
                    # XXX: There is some very strange behavior in tabix that causes it to
                    # return lines outside of the requested range. I can't figure out what
                    # causes it, but hopefully it has no practical effect. In this case,
                    # despite all depth matrix tables having md5-identical chrom and pos
                    # columns, (command line, so nothing to do with scan2's
                    # tabix functions) tabix <file> 1:2000001-3000000 returns
                    # position 2000000 for 2 out of the 17 matrix files (1278 and 5572) in
                    # in addition to the 950,000 other positions remaining after gap
                    # removal.
                    # This affects any tool using tabix as well as SCAN2. As a result, we
                    # just have to make sure the tilemap is restricted to the requested
                    # range and that these tabix data.tables are as well.
                    dpm.basepair <- read.tabix.data(f, region=chunks[i])
                    dpm.basepair <- dpm.basepair[pos >= start(chunks[i])[1] & pos <= end(chunks[i])[1],-(1:2)]
                    if (nrow(dpm.basepair) == 0) {
                        ret <- c()
                    } else {
                        # initialize to NA because all bp positions are not in the tile map
                        dpm.basepair[from(tilemap), tileid := to(tilemap)]
                        ret <- dpm.basepair[,setNames(as.list(colMeans(.SD)), colnames(.SD)),by=tileid][!is.na(tileid)][, -'tileid']
                    }
                })
                p(class='sticky', amount=1, pc)
                ret
            })
            mean.mat <- do.call(cbind, mean.mats)
        }
        list(tiles=tiles, tiles.not.in.gatk=tiles.not.in.gatk, mean.mat=mean.mat)
    })
}, enable=TRUE)

tiles <- do.call(c, lapply(results, function(r) r$tiles))
tiles.not.in.gatk <- do.call(c, lapply(results, function(r) r$tiles.not.in.gatk))
mean.mat <- rbindlist(lapply(results, function(r) r$mean.mat))
save(chrom, base.tile.size, tiles, tiles.not.in.gatk, mean.mat, file=out.rda)

cat('Final memory profile:\n')
print(gc())

if ('snakemake' %in% ls()) {
    sink()
}
