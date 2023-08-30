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
        snakemake@output['summary'],
        snakemake@input['qbeds']
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}


args <- commandArgs(trailingOnly=T)

if (length(args) < 5) {
    cat('NOTE: This script is only valid when applied to QBEDs with the same 1kb genome tiling as used in SCAN2\'s default spatial sensitivity model\n')
    stop("usage: enrichment_sensitivity_qbed.R scan2_full_object filtered_snvs.txt filtered_indels.txt summary_output.txt qbed_1kb1 [ qbed_1kb2 ... qbed_1kbN ]")
}

in.rda <- args[1]
in.filtered.snvs <- args[2]
in.filtered.indels <- args[3]
#out.full <- args[2]
out.summary <- args[4]
bed.files <- args[-(1:4)]

#for (f in c(out.full, out.summary)) {
for (f in c( out.summary)) {
    if (file.exists(f))
        stop(sprintf("output file %s already exists, please delete it first", f))
}

suppressMessages(library(scan2))
suppressMessages(library(mutenrich))


cat("loading SCAN2 object (must be a FULL object):", in.rda, "\n")
load(in.rda) # loads 'results'

# Frees ~6GB.
results@spatial.sensitivity$models <- NULL

# List of all mutations that survived the various duplicate and cluster filters
# Be sure to only keep the mutations for this sample.
filt.muts <- rbind(fread(in.filtered.snvs), fread(in.filtered.indels))[sample == results@single.cell]

called.muts <- results@gatk[pass == TRUE | rescue == TRUE]
called.muts[, retained := paste(chr,pos,refnt,altnt) %in% filt.muts[,paste(chr,pos,refnt,altnt)]]

# For restricting to autosomes
autosome.seqinfo <- seqinfo(genome.string.to.tiling(genome=results@genome.string, group='auto'))

# genome-wide tiling with 1kb resolution. Has all genomic covariates that we
# want summarized in the final output. Allows a simple join (or, really, even
# a cbind()) onto the QBED.
# IMPORTANT: as saved, sort is lexicographical!
# i.e., chr1, chr10, chr11, ... not chr1, chr2, chr3, ...
data <- results@spatial.sensitivity$data[!is.na(snv.n.training.neighborhood)]

# Make a GRanges version so we can get the limits of data regions
gdata <- GRanges(seqnames=data$chr, ranges=IRanges(start=data$start, end=data$end),
    seqinfo=autosome.seqinfo)

# IN OUR DATA, produces 1 range per chromosome, corresponds to non-NA values.
# doesn't produce 1 range per chromosome in general.
rgdata <- reduce(gdata)
if (length(rgdata) != length(autosome.seqinfo)) {
    print(as.data.frame(rgdata))
    stop('reduce(gdata) did not produce 22 contiguous ranges. see above.')
}


print(gc())

d.summary <- data.table()
for (i in 1:length(bed.files)) {
    bed.file <- bed.files[i]
    cat("processing QBED", i, "/", length(bed.files), bed.file, "\n")
    bed <- as(read.bed(bedfile=bed.file, genome=results@genome.seqinfo, is.qbed=TRUE), 'GRanges')
    metaline <- read.bed.metadata(bed.file, is.qbed=TRUE, return.line=TRUE)

    # Restrict BED file to just autosomes. Use this seqinfo object to enforce
    # uniformity across GRanges.
    seqlevels(bed, pruning.mode='coarse') <- seqlevels(autosome.seqinfo)
    seqinfo(bed) <- autosome.seqinfo

    # Restrict to regions present in `data` (which rgdata represents)
    bed <- restrict(bed,
        start=setNames(start(rgdata), seqnames(rgdata)),
        end=setNames(end(rgdata), seqnames(rgdata)))
    bed <- bed[width(bed) > 0,]  # somehow restrict() creates a few 0-width bins

    # Convert to bed (GRanges) to data.table so we can join
    dt <- data.table(sample=results@single.cell,
        chr=as.character(seqnames(bed)),
        start=start(bed),
        end=end(bed),
        width=width(bed),
        metaline=metaline,
        feature=bed$quantile)

    d.full <- dt[data,, on=.(chr,start,end)]
    # almost same ranges as the read.bed() result, but post-join so only those in
    # data and in the same order.
    bed2 <- GRanges(seqnames=d.full$chr, ranges=IRanges(start=d.full$start, end=d.full$end), seqinfo=autosome.seqinfo)
    for (mt in c('snv', 'indel')) {
        # use this for counting somatic calls
        cat(mt, 'somatic pass', '\n')
        sp <- count.somatic.sites.for.sens(grs=bed2,
            sites=called.muts[retained == TRUE & pass==TRUE & muttype==mt],
            seqinfo=seqinfo(bed2),
            neighborhood.tiles=10)[, .(n.calls)]
        colnames(sp) <- paste0('sum.', mt, '.pass.', colnames(sp))
    
        cat(mt, 'somatic rescue', '\n')
        sr <- count.somatic.sites.for.sens(grs=bed2,
            sites=called.muts[retained == TRUE & rescue==TRUE & muttype==mt],
            seqinfo=seqinfo(bed2),
            neighborhood.tiles=10)[, .(n.calls)]
        colnames(sr) <- paste0('sum.', mt, '.rescue.', colnames(sr))
    
        # We want it to be parallel, not just joinable
        if (any(d.full$chr != sp$chr))
            stop('some chr entries of d.full != sp')
        if (any(d.full$start != sp$start))
            stop('some start entries of d.full != sp')
        if (any(d.full$end != sp$end))
            stop('some end entries of d.full != sp')

        d.full <- cbind(d.full, sp, sr)
    }
    
    # Window weight per feature (sum(weight)=1 within each feature)
    d.full[, weight := width/sum(as.numeric(width)), by=.(metaline,feature)]
    this.summary <- d.full[,.(sample=sample[1],
        width=sum(as.numeric(width)),
        mean.abs.gp.mu=sum(weight*abs(gp.mu)),
        mean.gp.sd=sum(weight*gp.sd),
        mean.snv.n.training=sum(weight*snv.n.training),
        mean.snv.n.training.neighborhood=sum(weight*snv.n.training.neighborhood),
        sum.snv.n.training.passed=sum(snv.n.training.passed),
        sum.snv.n.training=sum(snv.n.training),
        snv.sens=sum(snv.n.training.passed)/sum(snv.n.training),
        sum.snv.n.pass=sum(sum.snv.pass.n.calls),
        sum.snv.n.rescue=sum(sum.snv.rescue.n.calls),
        sum.indel.n.training.passed=sum(indel.n.training.passed),
        sum.indel.n.training=sum(indel.n.training),
        indel.sens=sum(indel.n.training.passed)/sum(indel.n.training),
        sum.indel.n.pass=sum(sum.indel.pass.n.calls),
        sum.indel.n.rescue=sum(sum.indel.rescue.n.calls),
        mean.sc.dp=sum(weight*mean.sc.dp, na.rm=TRUE),
        mean.bulk.dp=sum(weight*mean.bulk.dp, na.rm=TRUE),
        sum.bases.gt.snv.sc.min.dp=sum(bases.gt.snv.sc.min.dp, na.rm=TRUE),
        sum.bases.gt.snv.bulk.min.dp=sum(bases.gt.snv.bulk.min.dp, na.rm=TRUE),
        sum.bases.gt.indel.sc.min.dp=sum(bases.gt.indel.sc.min.dp, na.rm=TRUE),
        sum.bases.gt.indel.bulk.min.dp=sum(bases.gt.indel.bulk.min.dp, na.rm=TRUE)),
        by=.(metaline,feature)][order(metaline,feature)]

    d.summary <- rbind(d.summary, this.summary)

    print(gc())
}
    
# No longer writing out the full file now that we're doing multiple QBEDs per run
#fwrite(d.full, file=out.full)
fwrite(d.summary, file=out.summary)

if ('snakemake' %in% ls()) {
    sink()
}
