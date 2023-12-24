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
        snakemake@input['scan2_tiles'],
        snakemake@input['tiles'],
        snakemake@output['csv']
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}


args <- commandArgs(trailingOnly=T)

if (length(args) != 6) {
    cat("input_tile_set must be a perfectly aligned cover of scan2_sensitivity_tiles.bed!\n")
    cat("the information in scan2_sensitivity_tiles.bed is not really used; it's already stored in the SCAN2 object. It is primarily here to force the user to acknowledge what SCAN2 tiles were used.\n")
    stop("usage: enrichment_sensitivity_for_tiles.R scan2_full_object filtered_snvs.txt filtered_indels.txt scan2_sensitivity_tiles.bed input_tile_set.bed output.csv")
}

in.rda <- args[1]
in.filtered.snvs <- args[2]
in.filtered.indels <- args[3]
scan2.tile.file <- args[4]
in.tile.file <- args[5]
out.csv <- args[6]

for (f in c(out.csv)) {
    if (file.exists(f))
        stop(sprintf("output file %s already exists, please delete it first", f))
}

suppressMessages(library(scan2))
suppressMessages(library(mutenrich))

cat("Reading tile files..\n")
gr2 <- function(tab) GenomicRanges::GRanges(seqnames=tab$chr, ranges=IRanges(start=tab$start, end=tab$end))
scan2.tiles <- fread(scan2.tile.file)
setnames(scan2.tiles, old=paste0("V",1:3), new=c('chr','start','end'))
tiles <- fread(in.tile.file)
setnames(tiles, old=paste0("V",1:3), new=c('chr','start','end'))
gtiles <- gr2(tiles)


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
scan2.tiles <- scan2.tiles[order(chr,start,end)] # matches the lexicographical sorting above
gscan2.tiles <- gr2(scan2.tiles)

# sanity check: assert that the 1kb tiles in `data` are indeed the 1kb tile set
if (any(scan2.tiles$chr != paste0('chr', data$chr))) {
    stop("'chr' columns differ between SCAN2 object and 1kb tiles")
} else if (any(scan2.tiles$start != data$start)) {
    stop("'start' columns differ between SCAN2 object and 1kb tiles")
} else if (any(scan2.tiles$end != data$end)) {
    stop("'end' columns differ between SCAN2 object and 1kb tiles")
} else {
    cat('sanity check PASSED: 1kb tile set matches tiles in SCAN2 object\n')
}

# this map connects the SCAN2 bins to the larger bins that are perfectly aligned
cat("Creating SCAN2 <-> input tile mapping\n")
ols <- findOverlaps(gscan2.tiles, gtiles)
print(str(ols))
#stop('done')
data[, kb.to.mb.map := to(ols)]
data[, w := end-start+1]  # for weighted means in next call

# Create the upscaled version of `data`. Sum count-type data or take a weighted
# mean for mean data.
data.new <- data[,.(sample=results@single.cell, chr=as.character(chr[1]), start=start[1], end=tail(end,1), width=tail(end,1)-start[1]+1, snv.n.training=sum(snv.n.training, na.rm=T), snv.n.training.passed=sum(snv.n.training.passed, na.rm=T), snv.n.training.neighborhood=sum(snv.n.training.neighborhood,na.rm=T), snv.n.calls=sum(snv.n.calls,na.rm=T), snv.n.neighborhood=sum(snv.n.neighborhood,na.rm=T), indel.n.training=sum(indel.n.training,na.rm=T), indel.n.training.passed=sum(indel.n.training.passed,na.rm=T), indel.n.training.neighborhood=sum(indel.n.training.neighborhood,na.rm=T), indel.n.calls=sum(indel.n.calls, na.rm=T), indel.n.neighborhood=sum(indel.n.neighborhood,na.rm=T), gp.mu=sum(gp.mu*w/sum(w,na.rm=T), na.rm=T), gp.sd=sum(gp.sd*w/sum(w,na.rm=T),na.rm=T), mean.sc.dp=sum(mean.sc.dp*w/sum(w,na.rm=T),na.rm=T), mean.bulk.dp=sum(mean.bulk.dp*w/sum(w,na.rm=T),na.rm=T), bases.gt.snv.sc.min.dp=sum(bases.gt.snv.sc.min.dp,na.rm=T), bases.gt.snv.bulk.min.dp=sum(bases.gt.snv.bulk.min.dp, na.rm=T), bases.gt.indel.sc.min.dp=sum(bases.gt.indel.sc.min.dp,na.rm=T), bases.gt.indel.bulk.min.dp=sum(bases.gt.indel.bulk.min.dp,na.rm=T)),
    by=kb.to.mb.map]

data <- data[!is.na(snv.n.training.neighborhood)]
data.new <- data.new[!is.na(snv.n.training.neighborhood)]  # does nothing for 1mb bins

gdata <- GRanges(seqnames=data$chr, ranges=IRanges(start=data$start, end=data$end),
    seqinfo=autosome.seqinfo)
gdata.new <- GRanges(seqnames=data.new$chr, ranges=IRanges(start=data.new$start, end=data.new$end),
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

    # repeat for larger bins
    g <- count.germline.sites.for.sens(grs=gdata.new,
        sites=results@gatk[training.site==TRUE & muttype==mt & dp >= dp.cutoff],
        seqinfo=seqinfo(gdata.new),
        neighborhood.tiles=10)[, .(n.training, n.training.passed)]  # this no longer makes sense, but is unused
    g[, sens := ifelse(n.training > 0, n.training.passed / n.training, 0)]
    colnames(g) <- paste0(mt, '.', colnames(g))
    # Update the counts in `data.new` and add the new <muttype>.sens column
    for (cn in colnames(g))
        data.new[[cn]] <- g[[cn]]
}

print(gc())

for (mt in c('snv', 'indel')) {
    # use this for counting somatic calls
    cat(mt, 'somatic pass', '\n')
    sp <- count.somatic.sites.for.sens(grs=gdata.new,
        sites=called.muts[retained == TRUE & pass==TRUE & muttype==mt],
        seqinfo=seqinfo(gdata.new),
        neighborhood.tiles=10)[, .(n.calls)]  # again, doesn't make sense but is unused
    colnames(sp) <- paste0('sum.', mt, '.pass.', colnames(sp))

    cat(mt, 'somatic rescue', '\n')
    sr <- count.somatic.sites.for.sens(grs=gdata.new,
        sites=called.muts[retained == TRUE & rescue==TRUE & muttype==mt],
        seqinfo=seqinfo(gdata.new),
        neighborhood.tiles=10)[, .(n.calls)]
    colnames(sr) <- paste0('sum.', mt, '.rescue.', colnames(sr))

    # We want it to be parallel, not just joinable
    if (any(data.new$chr != sp$chr))
        stop('some chr entries of data.new != sp')
    if (any(data.new$start != sp$start))
        stop('some start entries of data.new != sp')
    if (any(data.new$end != sp$end))
        stop('some end entries of data.new != sp')

    data.new <- cbind(data.new, sp, sr)
    
    print(gc())
}
    
fwrite(data.new, file=out.csv)

if ('snakemake' %in% ls()) {
    sink()
}
