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
        snakemake@output['summary'],
        snakemake@input['beds']
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}


args <- commandArgs(trailingOnly=T)

if (length(args) < 3) {
    cat("WARNING: this script manually removes all features on non-autosomes (chrs1-22)\n")
    stop("usage: enrichment_sensitivity_qbed.R scan2_full_object summary_output bed1 [ bed2 ... bed N ]")
}

in.rda <- args[1]
out.summary <- args[2]
bed.files <- args[-(1:2)]

for (f in c(out.summary)) {
    if (file.exists(f))
        stop(sprintf("output file %s already exists, please delete it first", f))
}

suppressMessages(library(scan2))
suppressMessages(library(mutenrich))

cat("loading SCAN2 object (must be a FULL object):", in.rda, "\n")
load(in.rda) # loads 'results'

d.summary <- data.table()
for (i in 1:length(bed.files)) {
    bed.file <- bed.files[i]
    cat("processing BED", i, "/", length(bed.files), bed.file, "\n")
    bed <- read.bed(bedfile=bed.file, genome=results@genome.seqinfo, is.qbed=FALSE, has.metaline=TRUE)
    metaline <- read.bed.metadata(bed.file, is.qbed=FALSE, return.line=TRUE)

    # only keep autosomal features
    bed <- bed[seqnames(bed) %in% 1:22,]

    s <- countOverlaps(bed, gr(results@gatk[muttype == 'snv' & training.site == TRUE]))
    sp <- countOverlaps(bed, gr(results@gatk[muttype == 'snv' & training.site == TRUE & training.pass == TRUE]))
    
    i <- countOverlaps(bed, gr(results@gatk[muttype == 'indel' & training.site == TRUE]))
    ip <- countOverlaps(bed, gr(results@gatk[muttype == 'indel' & training.site == TRUE & training.pass == TRUE]))
    
    d.full <- rbind(data.table(
            chr=as.character(seqnames(bed)),
            start=start(bed),
            end=end(bed),
            muttype='snv',
            width=width(bed),
            quantile=bed$feature,  # call it quantile for compatibility with qbeds
            total=s,
            pass=sp),
        data.table(
            chr=as.character(seqnames(bed)),
            start=start(bed),
            end=end(bed),
            muttype='indel',
            width=width(bed),
            quantile=bed$feature,
            total=i,
            pass=ip)
    )
    
    this.summary <- d.full[,.(meta=metaline, sample=results@single.cell, size.mb=sum(width)/1e6, pass=sum(pass), total=sum(total), sens=sum(pass)/sum(total)),by=.(muttype,quantile)][order(muttype,quantile)]

    d.summary <- rbind(d.summary, this.summary)
}
    
# No longer writing out the full file now that we're doing multiple QBEDs per run
#fwrite(d.full, file=out.full)
fwrite(d.summary, file=out.summary)

if ('snakemake' %in% ls()) {
    sink()
}
