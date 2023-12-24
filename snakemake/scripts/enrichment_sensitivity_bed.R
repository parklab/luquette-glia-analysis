#!/usr/bin/env Rscript

# detect script being run by snakemake
# if so, make a mock commandArgs function
if ('snakemake' %in% ls()) {
    logfile <- snakemake@log[[1]]
    con <- file(logfile, 'w')
    sink(con, type='output')
    sink(con, type='message')

    commandArgs <- function(...) unlist(c(
        snakemake@input['scan2_full_object'],
        snakemake@input['filtered_snvs'],
        snakemake@input['filtered_indels'],
        snakemake@output['full'],
        snakemake@output['summary'],
        snakemake@input[['beds']]
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}


args <- commandArgs(trailingOnly=T)

if (length(args) < 6) {
    cat("WARNING: this script manually removes all features on non-autosomes (chrs1-22)\n")
    stop("usage: enrichment_sensitivity_bed2.R scan2_full_object filtered_snvs.txt filtered_indels.txt full_output.txt summary_output.txt bed1 [ bed2 ... bedN ]")
}

in.rda <- args[1]
in.filtered.snvs <- args[2]
in.filtered.indels <- args[3]
out.full <- args[4]
out.summary <- args[5]
bed.files <- args[-(1:5)]

for (f in c(out.full, out.summary)) {
    if (file.exists(f))
        stop(sprintf("output file %s already exists, please delete it first", f))
}

suppressMessages(library(scan2))
suppressMessages(library(mutenrich))

cat("loading SCAN2 object (must be a FULL object):", in.rda, "\n")
load(in.rda) # loads 'results'

# ~6GB of models. doesn't reduce peak usage, but helps to ensure that metadata
# calculation doesn't exceed memlimit.
results@spatial.sensitivity$models <- NULL

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

# List of all mutations that survived the various duplicate and cluster filters
# Be sure to only keep the mutations for this sample.  This filtration is not
# reflected in the SCAN2 object because it's performed in a post-processing script
# that aggregates data from all single cells.
filt.muts <- rbind(fread(in.filtered.snvs), fread(in.filtered.indels))[sample == results@single.cell]

called.muts <- results@gatk[pass == TRUE | rescue == TRUE]
called.muts[, retained := paste(chr,pos,refnt,altnt) %in% filt.muts[,paste(chr,pos,refnt,altnt)]]


# Restrict BED file to just autosomes. Use this seqinfo object to enforce
# uniformity across GRanges.
autosome.seqinfo <- seqinfo(genome.string.to.tiling(genome=results@genome.string, group='auto'))
seqlevels(gbed, pruning.mode='coarse') <- seqlevels(autosome.seqinfo)
seqinfo(gbed) <- autosome.seqinfo
gbed <- sort(gbed)
print(gc())


# For AB model parameters (gp.mu, gp.sd, number of training sites in the local
# neighborhood), use the 1kb resolution estimates already made as part of the
# general spatial model.
#
# These should NOT be recomputed for each window because for, e.g., large windows
# like introns the midpoint AB calculation will not be representative of the AB
# across the window. Based on the length of local AB correlation in PTA products,
# the 1kb windows used for the general spatial sensitivity model ARE reasonable
# representations of the whole window. Speculatively, the largest violation of
# this assumption is probably when a phasing switch error occurs in an imbalanced
# region and causes the causes an abrupt switch in local AB estimate.
ss <- results@spatial.sensitivity$data
gdata <- GRanges(seqnames=ss$chr,
    ranges=IRanges(start=ss$start, end=ss$end),
    seqinfo=autosome.seqinfo)
gdata$abs.gp.mu <- abs(ss$gp.mu)
gdata$gp.sd <- ss$gp.sd
# N.B. although indel training sites exist, they are not used for estimating AB,
# even for indel calling.
gdata$snv.n.training <- ss$snv.n.training
gdata$snv.n.training.neighborhood <- ss$snv.n.training.neighborhood
gdata <- sort(gdata)

gdata <- gdata[!is.na(gdata$snv.n.training.neighborhood),]

# IN OUR DATA, produces 1 range per chromosome, corresponds to non-NA values.
# doesn't produce 1 range per chromosome in general.
rgdata <- reduce(gdata)
if (length(rgdata) != length(autosome.seqinfo)) {
    print(as.data.frame(rgdata))
    stop('reduce(gdata) did not produce 22 contiguous ranges. see above.')
}

gbed <- restrict(gbed,
    start=setNames(start(rgdata), seqnames(rgdata)),
    end=setNames(end(rgdata), seqnames(rgdata)))


# Append binned averages as metadata to `gbed`
for (signal.name in c('abs.gp.mu', 'gp.sd', 'snv.n.training', 'snv.n.training.neighborhood')) {
#for (signal.name in c('abs.gp.mu')) {
    cat(paste0("computing binned averages for signal=", signal.name, '\n'))
    print(system.time(rlelist <- mcolAsRleList(gdata, signal.name)))
    cat('binnedAverage\n')
    # XXX: na.rm=TRUE increases runtime from 0.115 seconds -> 1215.156 seconds.
    # nope. that's not a typo. it's really a 10,000-fold increase in runtime for 19
    # ranges with NA values out of 572784 ranges in a test dataset.  each of those 19
    # windows is 1kb except one window of 566bp.
    print(system.time(ba <- binnedAverage(gbed, rlelist, signal.name))) #, na.rm=TRUE)))
    cat('setting mcols()\n')
    print(system.time(mcols(gbed)[[signal.name]] <- mcols(ba)[[signal.name]]))
}

cat("Counting germline and somatic calls per feature\n")
for (mt in c('snv', 'indel')) {
    cat(mt, 'germline', '\n')
    # XXX: it would be ideal to also use the bulk cutoff, since this is used for permuting.
    # however, our depth covariate analysis currently only collects %bases >= sc.min.dp and
    # %bases >= bulk.min.dp separately, not %bases >= sc.min.dp & bulk.min.dp. Approximating
    # by sc.min.dp should be very close to correct since it is rare for bulk.dp < sc.dp
    dp.cutoff <- results@static.filter.params[[mt]]$min.sc.dp
    g <- count.germline.sites.for.sens(grs=gbed,
        sites=results@gatk[training.site==TRUE & muttype==mt & dp >= dp.cutoff],
        seqinfo=seqinfo(gbed),
        neighborhood.tiles=10)[, .(n.training, n.training.passed)]
    colnames(g) <- paste0('sum.', mt, '.', colnames(g))

    cat(mt, 'somatic pass', '\n')
    sp <- count.somatic.sites.for.sens(grs=gbed,
        sites=called.muts[retained == TRUE & pass==TRUE & muttype==mt],
        seqinfo=seqinfo(gbed),
        neighborhood.tiles=10)[, .(n.calls)]
    colnames(sp) <- paste0('sum.', mt, '.pass.', colnames(sp))

    cat(mt, 'somatic rescue', '\n')
    sr <- count.somatic.sites.for.sens(grs=gbed,
        sites=called.muts[retained == TRUE & rescue==TRUE & muttype==mt],
        seqinfo=seqinfo(gbed),
        neighborhood.tiles=10)[, .(n.calls)]
    colnames(sr) <- paste0('sum.', mt, '.rescue.', colnames(sr))

    mcols(gbed) <- cbind(mcols(gbed), g, sp, sr)
}

cat("Final GRanges object with binned signals:\n")
print(gbed)

d.full <- data.table(sample=results@single.cell, width=width(gbed), as.data.table(mcols(gbed)))
# Window weight per feature (sum(weight)=1 within each feature)
d.full[, weight := width/sum(as.numeric(width)), by=.(metaline,feature)]
d.summary <- d.full[,.(sample=sample[1],
        width=sum(as.numeric(width)),
        mean.abs.gp.mu=sum(weight*abs.gp.mu),
        mean.gp.sd=sum(weight*gp.sd),
        mean.snv.n.training=sum(weight*snv.n.training),
        mean.snv.n.training.neighborhood=sum(weight*snv.n.training.neighborhood),
        sum.snv.n.training.passed=sum(sum.snv.n.training.passed),
        sum.snv.n.training=sum(sum.snv.n.training),
        snv.sens=sum(sum.snv.n.training.passed)/sum(sum.snv.n.training),
        sum.snv.n.pass=sum(sum.snv.pass.n.calls),
        sum.snv.n.rescue=sum(sum.snv.rescue.n.calls),
        sum.indel.n.training.passed=sum(sum.indel.n.training.passed),
        sum.indel.n.training=sum(sum.indel.n.training),
        indel.sens=sum(sum.indel.n.training.passed)/sum(sum.indel.n.training),
        sum.indel.n.pass=sum(sum.indel.pass.n.calls),
        sum.indel.n.rescue=sum(sum.indel.rescue.n.calls)
    ), by=.(metaline,feature)][order(metaline,feature)]

# No longer writing out the full file now that we're doing multiple QBEDs per run
fwrite(d.full, file=out.full)
fwrite(d.summary, file=out.summary)

print(gc())

if ('snakemake' %in% ls()) {
    sink()
}
