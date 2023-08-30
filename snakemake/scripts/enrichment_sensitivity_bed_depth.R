#!/usr/bin/env Rscript

# detect script being run by snakemake
# if so, make a mock commandArgs function
if ('snakemake' %in% ls()) {
    logfile <- snakemake@log[[1]]
    con <- file(logfile, 'w')
    sink(con, type='output')
    sink(con, type='message')

    commandArgs <- function(...) unlist(c(
        #snakemake@params['n_cores'],  # don't pass this as a param or else it triggers reruns on "param change"
        as.integer(snakemake@threads),
        snakemake@params['sample'],
        snakemake@params['bulk'],
        snakemake@params['genome'],
        snakemake@input['yaml'],
        snakemake@input['dptab'],
        snakemake@output['full'],
        snakemake@output['summary'],
        snakemake@input[['beds']]
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}


args <- commandArgs(trailingOnly=T)

if (length(args) < 9) {
    cat("Do not use this script to calculate depth over large numbers of fine tilings (like 1000bp qbeds or 200bp chromHMM annotations). The memory usage will be too high to multithread effectively.\n")
    cat("WARNING: this script manually removes all features on non-autosomes (chrs1-22)\n")
    stop("usage: enrichment_sensitivity_bed2.R n_cores single_cell_sample_ID bulk_ID genome_string scan2_config.yaml joint_depth_table.tab.gz full_output.txt summary_output.txt bed1 [ bed2 ... bedN ]")
}

n.cores <- as.integer(args[1])
single.cell.id <- args[2]
bulk.id <- args[3]
genome.string <- args[4]
config.file <- args[5]
dptab.file <- args[6]
out.full <- args[7]
out.summary <- args[8]
bed.files <- args[-(1:8)]

#testing
#n.cores <- 16
#single.cell.id <- "UMB4976_E1"
#bulk.id <- "4976-190613-cer"
#genome.string <- 'hs37d5'
#config.file <- 'scan2/UMB4976/scan2/scan.yaml'
#dptab.file <- 'scan2/UMB4976/scan2/depth_profile/joint_depth_matrix.tab.gz'
#bed.files <- c('enrichment/nott/bed_regions/masked_bed/nott_enhprom___overlapping_peaks.bed',
#'enrichment/gencode_simplified/bed_regions/masked_bed/gencode_simplified___overlapping_peaks.bed',
#'enrichment/gtex_genes/bed_regions/masked_bed/gtex___genes.bed',
#'enrichment/roadmap_chromhmm_brain/bed_regions/masked_bed/roadmap___chromhmm___15state___E073.bed')

for (f in c(out.full, out.summary)) {
    if (file.exists(f))
        stop(sprintf("output file %s already exists, please delete it first", f))
}

suppressMessages(library(scan2))
suppressMessages(library(mutenrich))

if (n.cores > 1) {
    cat("Using library(future) with", n.cores, "cores\n")
    suppressMessages(library(future))
    plan(multicore, workers=n.cores)
}

# If we read in the actual SCAN2 object then a large amount of memory will be irreversibly
# allocated, making future-style (fork) multicore infeasible.
# So instead, build up a minimal SCAN2 object with no internal data.
results <- make.scan(single.cell=single.cell.id, bulk=bulk.id, genome=genome.string)
results <- add.static.filter.params(results, config.path=config.file)


# Create a single GRanges file containing all BED intervals.
# Prepend features from each BED file with that BED's metadata line.
gbed <- GRanges()
for (i in 1:length(bed.files)) {
    bed.file <- bed.files[i]
    cat("processing BED", i, "/", length(bed.files), bed.file, "\n")
    new.gbed <- as(read.bed(bed.file, is.qbed=FALSE, has.metaline=TRUE, genome=results@genome.seqinfo), 'GRanges')
    new.gbed$metaline <- read.bed.metadata(bed.file, is.qbed=FALSE, return.line=TRUE)
    gbed <- c(gbed, new.gbed)
    print(gc())
}

# Restrict BED file to just autosomes. Use this seqinfo object to enforce
# uniformity across GRanges.
autosome.seqinfo <- seqinfo(genome.string.to.tiling(genome=results@genome.string, group='auto'))
seqlevels(gbed, pruning.mode='coarse') <- seqlevels(autosome.seqinfo)
seqinfo(gbed) <- autosome.seqinfo
gbed <- sort(gbed)


# For depth covariates, compute averages by directly accessing basepair-level data.
# The 1MB parallelization tiling uses ~2.5GB at peak per thread.
#
# I don't know why
# it doesn't scale with window size, but 2.5GB is peak usage for 100kB and only hits
# ~3.5GB for 5MB tilings.
# This memory non-scaling even happens with an empty SCAN2 object (containing just
# single cell ID, bulk ID, static.filter.params) as created above.
#
# If a 200bp resolution BED is added, peak usage is ~4G/thread.
#options(future.globals.maxSize=5 * 1024^3)
d <- compute.spatial.sensitivity.depth(
    single.cell.id=results@single.cell,
    bulk.id=results@bulk,
    static.filter.params=results@static.filter.params,
    joint.dptab.path=dptab.file,
    genome.string=results@genome.string,
    grs.for.sens=gbed,
    quiet=FALSE,
    report.mem=FALSE)


d.full <- data.table(sample=results@single.cell, width=width(gbed), as.data.table(mcols(gbed)), d[,-(1:3)])[!is.na(mean.sc.dp)]
# Window weight per feature (sum(weight)=1 within each feature)
d.full[, weight := width/sum(as.numeric(width)), by=.(metaline,feature)]
d.summary <- d.full[,.(sample=sample[1],
        width=sum(as.numeric(width)),
        mean.sc.dp=sum(weight*mean.sc.dp),
        mean.bulk.dp=sum(weight*mean.bulk.dp),
        sum.bases.gt.snv.sc.min.dp=sum(bases.gt.snv.sc.min.dp),
        sum.bases.gt.snv.bulk.min.dp=sum(bases.gt.snv.bulk.min.dp),
        sum.bases.gt.indel.sc.min.dp=sum(bases.gt.indel.sc.min.dp),
        sum.bases.gt.indel.bulk.min.dp=sum(bases.gt.indel.bulk.min.dp)
    ), by=.(metaline,feature)][order(metaline,feature)]

# No longer writing out the full file now that we're doing multiple QBEDs per run
fwrite(d.full, file=out.full)
fwrite(d.summary, file=out.summary)

if ('snakemake' %in% ls()) {
    sink()
}
