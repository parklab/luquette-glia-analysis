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
        snakemake@input['tiles_1kb'],
        snakemake@input['tiles_1mb'],
        snakemake@output['summary'],
        snakemake@input['qbeds']
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}


args <- commandArgs(trailingOnly=T)

if (length(args) < 7) {
    cat('NOTE: This script is only valid when applied to QBEDs with the same 1kb genome tiling as used in SCAN2\'s default spatial sensitivity model OR QBEDs with a specific 1MB genome tiling that aligns to the 1kb tiling from the sensitivity model\n')
    stop("usage: enrichment_sensitivity_qbed_1kb_1mb.R scan2_full_object filtered_snvs.txt filtered_indels.txt tile_set_1kb.bed tile_set_1mb.bed summary_output.txt qbed1 [ qbed2 ... qbedN ]")
}

in.rda <- args[1]
in.filtered.snvs <- args[2]
in.filtered.indels <- args[3]
in.tile.1kb.file <- args[4]
in.tile.1mb.file <- args[5]
#out.full <- args[2]
out.summary <- args[6]
bed.files <- args[-(1:6)]

#for (f in c(out.full, out.summary)) {
for (f in c( out.summary)) {
    if (file.exists(f))
        stop(sprintf("output file %s already exists, please delete it first", f))
}

suppressMessages(library(scan2))
suppressMessages(library(mutenrich))

cat("Reading tile files..\n")
gr2 <- function(tab) GenomicRanges::GRanges(seqnames=tab$chr, ranges=IRanges(start=tab$start, end=tab$end))
tiles.1kb <- fread(in.tile.1kb.file)
setnames(tiles.1kb, old=paste0("V",1:3), new=c('chr','start','end'))
tiles.1mb <- fread(in.tile.1mb.file)
setnames(tiles.1mb, old=paste0("V",1:3), new=c('chr','start','end'))
gtiles.1mb <- gr2(tiles.1mb)


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
data <- results@spatial.sensitivity$data
tiles.1kb <- tiles.1kb[order(chr,start,end)] # matches the lexicographical sorting above
gtiles.1kb <- gr2(tiles.1kb)

# sanity check: assert that the 1kb tiles in `data` are indeed the 1kb tile set
if (any(tiles.1kb$chr != paste0('chr', data$chr))) {
    stop("'chr' columns differ between SCAN2 object and 1kb tiles")
} else if (any(tiles.1kb$start != data$start)) {
    stop("'start' columns differ between SCAN2 object and 1kb tiles")
} else if (any(tiles.1kb$end != data$end)) {
    stop("'end' columns differ between SCAN2 object and 1kb tiles")
} else {
    cat('sanity check PASSED: 1kb tile set matches tiles in SCAN2 object\n')
}

# this map connects the 1kb bins to the larger 1mb bins that are perfectly aligned
cat("Creating 1kb <-> 1mb tile mapping\n")
ols <- findOverlaps(gtiles.1kb, gtiles.1mb)
print(str(ols))
#stop('done')
data[, kb.to.mb.map := to(ols)]
data[, w := end-start+1]  # for weighted means in next call

# Create the 1mb version of `data`. Sum count-type data or take a weighted
# mean for mean data.
data.1mb <- data[,.(chr=chr[1], start=start[1], end=tail(end,1), snv.n.training=sum(snv.n.training, na.rm=T), snv.n.training.passed=sum(snv.n.training.passed, na.rm=T), snv.n.training.neighborhood=sum(snv.n.training.neighborhood,na.rm=T), snv.n.calls=sum(snv.n.calls,na.rm=T), snv.n.neighborhood=sum(snv.n.neighborhood,na.rm=T), indel.n.training=sum(indel.n.training,na.rm=T), indel.n.training.passed=sum(indel.n.training.passed,na.rm=T), indel.n.training.neighborhood=sum(indel.n.training.neighborhood,na.rm=T), indel.n.calls=sum(indel.n.calls, na.rm=T), indel.n.neighborhood=sum(indel.n.neighborhood,na.rm=T), gp.mu=sum(gp.mu*w/sum(w,na.rm=T), na.rm=T), gp.sd=sum(gp.sd*w/sum(w,na.rm=T),na.rm=T), mean.sc.dp=sum(mean.sc.dp*w/sum(w,na.rm=T),na.rm=T), mean.bulk.dp=sum(mean.bulk.dp*w/sum(w,na.rm=T),na.rm=T), bases.gt.snv.sc.min.dp=sum(bases.gt.snv.sc.min.dp,na.rm=T), bases.gt.snv.bulk.min.dp=sum(bases.gt.snv.bulk.min.dp, na.rm=T), bases.gt.indel.sc.min.dp=sum(bases.gt.indel.sc.min.dp,na.rm=T), bases.gt.indel.bulk.min.dp=sum(bases.gt.indel.bulk.min.dp,na.rm=T)),
    by=kb.to.mb.map]

data <- data[!is.na(snv.n.training.neighborhood)]
data.1mb <- data.1mb[!is.na(snv.n.training.neighborhood)]  # does nothing for 1mb bins

gdata <- GRanges(seqnames=data$chr, ranges=IRanges(start=data$start, end=data$end),
    seqinfo=autosome.seqinfo)
gdata.1mb <- GRanges(seqnames=data.1mb$chr, ranges=IRanges(start=data.1mb$start, end=data.1mb$end),
    seqinfo=autosome.seqinfo)

# Need to re-score germline het SNP sensitivity only on sites with sufficient DP.
# This allows differentiation of sensitivity due to lack of DP with other reasons.
# In particular, this is important for correcting our enrichment analyses for
# sensitivity since some regions are very skewed for low DP (e.g., the "excluded"
# region is often ~33% > min.dp while others are >=80-90%). Using all hSNPs mixes
# the effect of low DP and 0 DP sensitivity; the majority of the latter case (0 DP)
# is in unassembled genomic regions where there are *ALSO* no hSNPs to score.
for (mt in c('snv', 'indel')) {
    dp.cutoff <- results@static.filter.params[[mt]]$min.sc.dp
    g <- count.germline.sites.for.sens(grs=gdata,
        # XXX: As is the case for BED file sensitivity analysis, it would be ideal to
        # require dp >= min.sc.dp & bulk.dp >= min.bulk.dp, but our %bases >= cutoff
        # only does each of these separately, not jointly.  It'd be good to fix this
        # in the future, but even if we don't, the approximation should be very close
        # since, in unassembled genomic regions, sc.dp=bulk.dp=0 and, in other regions,
        # sc.dp < bulk.dp is usually the case.
        sites=results@gatk[training.site==TRUE & muttype==mt & dp >= dp.cutoff],
        seqinfo=seqinfo(gdata),
        neighborhood.tiles=10)[, .(n.training, n.training.passed)]
    colnames(g) <- paste0(mt, '.', colnames(g))
    # Update the counts in `data`
    for (cn in colnames(g))
        data[[cn]] <- g[[cn]]

    # repeat for 1mb bins
    g <- count.germline.sites.for.sens(grs=gdata.1mb,
        sites=results@gatk[training.site==TRUE & muttype==mt & dp >= dp.cutoff],
        seqinfo=seqinfo(gdata.1mb),
        neighborhood.tiles=10)[, .(n.training, n.training.passed)]
    colnames(g) <- paste0(mt, '.', colnames(g))
    # Update the counts in `data`
    for (cn in colnames(g))
        data.1mb[[cn]] <- g[[cn]]
}


# Get the limits of data regions.
# This is NOT applied to 1mb bins.
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
    # the metaline contains all QBED info such as BINSIZE and number of quantiles.
    # those can be parsed out later as needed. packing these into one column allows
    # QBED sensitivity tables to be combined with regular BED tables.
    metaline <- read.bed.metadata(bed.file, is.qbed=TRUE, return.line=TRUE)
    # this call parses the metadata line so the BINSIZE element can be extracted
    metadata <- read.bed.metadata(bed.file, is.qbed=TRUE, return.line=FALSE)

    # Restrict BED file to just autosomes. Use this seqinfo object to enforce
    # uniformity across GRanges.
    seqlevels(bed, pruning.mode='coarse') <- seqlevels(autosome.seqinfo)
    seqinfo(bed) <- autosome.seqinfo

    if (metadata$BINSIZE == '1000') {
        cat(paste0('detected 1kb binning for bed file=', bed.file, '\n'))
        to.join <- data
        # Restrict to regions present in `data` (which rgdata represents)
        bed <- restrict(bed,
            start=setNames(start(rgdata), seqnames(rgdata)),
            end=setNames(end(rgdata), seqnames(rgdata)))
        bed <- bed[width(bed) > 0,]  # somehow restrict() creates a few 0-width bins
    } else if (metadata$BINSIZE == '1000000') {
        cat(paste0('detected 1Mb binning for bed file=', bed.file, '\n'))
        to.join <- data.1mb
    } else {
        stop(paste0('got BINSIZE=', metadata$BINSIZE, ' from QBED file=', bed.file, '. Expected either 1000 or 1000000.'))
    }

    # Convert to bed (GRanges) to data.table so we can join
    dt <- data.table(sample=results@single.cell,
        chr=as.character(seqnames(bed)),
        start=start(bed),
        end=end(bed),
        width=width(bed),
        metaline=metaline,
        feature=bed$quantile)
    d.full <- dt[to.join,,on=.(chr,start,end)]

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
