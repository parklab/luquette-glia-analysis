#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)

if (length(args) != 2)
    stop('usage: make_gene_regions.R gencode.v26lift37.annotation.gtf out.bed')


gencode.in <- args[1]
out.bed <- args[2]

if (file.exists(out.bed))
    stop(paste('output file', out.bed, 'already exists, please delete it first'))


library(rtracklayer)
library(GenomicRanges)
library(GenomeInfoDb)
library(data.table)

gn <- import.gff(gencode.in)

# Priority:
#   1. cds - Coding sequence
#   2. utr3 and utr5 - UTRs
#   3. intron - Non-coding, non-UTR regions in protein coding genes (labelled introns)
#   4. other - Everything else in the gene model (not sure what to call these)
#   5. upstream/downstream - Defer to other transcripts first, but otherwise differentiate
#          immediate up/downstream from other intergenic regions. Currently using 1 kb.
#   5. intergenic - Everything not in the gene model or gene proximal.

cds <- reduce(gn[gn$type=='CDS',], ignore.strand=T)
# > sum(width(cds))/1e6
# [1] 35.53679


utrs <- gn[gn$type=='UTR',]
txs <- gn[gn$type=='transcript',]

# Get a fast mapping between UTR transcript_id and transcript transcript_id
tmap <- data.table(idx=1:length(txs), transcript_id=txs$transcript_id)
setkey(tmap, transcript_id)
umap <- tmap[utrs$transcript_id, idx]

# Assign 5' and 3' UTR status
# A majority of sequence space has either left=0 or right=0, but
# surprisingly, almost 40% of all transcripts do not map the UTR
# to the first or last bases.
left <- start(utrs) - start(txs)[umap] # positive distance from tx start to UTR start
right <- end(txs)[umap] - end(utrs)    # positive distance from tx end to UTR end

utrs$utr.type <- ifelse(left <= right & strand(utrs) == '+' |
                   right <= left & strand(utrs) == '-', '5_UTR', '3_UTR')

# This isn't perfect, but the stats suggest it is doing a good job.
# The average size of a 5' UTR is 137 bp (vs. reported 100-200 for humans)
# and the average size of a 3' UTR is 508 bp (vs. reported ~800).
# > sapply(split(width(utrs), utrs$utr.type), sum)/1e6
#    3_UTR    5_UTR 
# 75.01655 18.83736    # NOTE: non-reduced
# > sapply(split(width(utrs), utrs$utr.type), mean)
#    3_UTR    5_UTR 
# 508.6212 137.7131 

utrs3 <- utrs[utrs$utr.type=='3_UTR',]
utrs5 <- utrs[utrs$utr.type=='5_UTR',]

# There's a small (<1MB overlap between the two sets). I'm just going to
# remove it since there doesn't seem to be a better way.
# > sum(width(reduce(utrs3)))/1e6
# [1] 43.16909      # NOTE: did not include ignore.strand=T in either of these
# > sum(width(reduce(utrs5)))/1e6
# [1] 11.08784
# > sum(width(intersect(reduce(utrs5), reduce(utrs3))))/1e6
# [1] 0.830637
utrs3 <- reduce(utrs3, ignore.strand=T)
utrs5 <- reduce(utrs5, ignore.strand=T)
unclassified <- intersect(utrs3, utrs5)

utrs3 <- setdiff(utrs3, unclassified)
utrs5 <- setdiff(utrs5, unclassified)
# > sum(width(utrs5))/1e6
# [1] 10.12833
# > sum(width(utrs3))/1e6
# [1] 41.85798

# Loci with any CDS annotation are preferred over UTR annotations
utrs3 <- setdiff(utrs3, cds)
utrs5 <- setdiff(utrs5, cds)
# > sum(width(utrs3))/1e6
# [1] 36.40233
# > sum(width(utrs5))/1e6
# [1] 8.601558

# Keep running total of what has already been claimed
covered <- union(cds, union(utrs3, utrs5, ignore.strand=T), ignore.strand=T)
# > sum(width(covered))/1e6
# [1] 80.54068

introns <- reduce(gn[gn$type == 'gene' & gn$gene_type == 'protein_coding',], ignore.strand=T)
# As expected, currently covered regions are 100% overlap with
# 'introns' so far, since no attempt has been made at making introns.
# > sum(width(introns))/1e6
# [1] 1276.634
# > sum(width(setdiff(introns,covered)))/1e6
# [1] 1196.286
introns <- setdiff(introns, covered, ignore.strand=T)

covered <- union(covered, introns, ignore.strand=T)
# > sum(width(covered))/1e6
# [1] 1276.827

# Everything else
other <- setdiff(reduce(gn, ignore.strand=T), covered, ignore.strand=T)
# > sum(width(other))/1e6
# [1] 326.4511
covered <- union(covered, other, ignore.strand=T)


# Generate 1kb up and downstream regions from each gene.
window <- 1000
p <- gn[strand(gn) == '+',]
n <- gn[strand(gn) == '-',]
up <- GRanges(seqnames=c(seqnames(p), seqnames(n)),
    ranges=IRanges(start=c(start(p) - window, end(n) + 1),
        end=c(start(p)-1, end(n) + window)))
up <- reduce(up, ignore.strand=T)
up <- setdiff(up, covered, ignore.strand=T)
down <- GRanges(seqnames=c(seqnames(p), seqnames(n)),
    ranges=IRanges(start=c(end(p) + 1, start(n) - window),
        end=c(end(p)+window, start(n) - 1)))
down <- reduce(down, ignore.strand=T)
down <- setdiff(down, covered, ignore.strand=T)
# After removing previously covered regions, the up and downstream
# sets don't overlap much. Just arbitrarily assign this small set to
# upstream.
# > sum(width(up))/1e6
# [1] 30.63927
# > sum(width(down))/1e6
# [1] 32.54219
# > sum(width(intersect(up,down,ignore.strand=T)))/1e6
# [1] 1.888588
down <- setdiff(down, up)

covered <- union(covered, union(up, down, ignore.strand=T), ignore.strand=T)

# Make final ranges with classifications
cds$type <- 'cds'
utrs3$type <- 'utr3'
utrs5$type <- 'utr5'
introns$type <- 'intron'
other$type <- 'other'
up$type <- 'upstream'
down$type <- 'downstream'
final <-  c(cds, utrs3, utrs5, introns, other, up, down)
# > sapply(split(width(final), final$type),sum)/1e6
#         cds  downstream      intron       other    upstream        utr3 
#   35.536790   30.653603 1196.286354  326.451135   30.639266   36.402334 
#        utr5 
#    8.601558 
final <- final[seqnames(final) %in% paste0('chr',c(1:22,'X','Y')),]


# Get intergenic spaces. Properly doing this requires getting each
# chromosome legnth. The correct genome to use is GRCh37, but because
# we use the 'chr' prefixes, hg19 is more convenient. The primary
# assembly sizes are nearly identical.
seqlevels(final) <- seqlevels(GenomeInfoDb::Seqinfo(genome='hg19'))  # yes, this has to happen before seqinfo<-
seqinfo(final) <- GenomeInfoDb::Seqinfo(genome='hg19')

# gaps() treats +/-/* strands all separately
intergenic <- gaps(final)
intergenic <- intergenic[strand(intergenic)=='*',]
intergenic <- intergenic[seqnames(intergenic) %in% paste0('chr',c(1:22,'X','Y')),]
intergenic$type <- 'intergenic'
# > sum(width(intergenic))/1e6
# [1] 1431.487

final <- c(final, intergenic)
final <- sort(final)
# > sum(width(final))/1e6
# [1] 3095.677
# > sapply(split(width(final), final$type),sum)/1e6
#         cds  downstream  intergenic      intron       other    upstream 
#   35.525457   30.625055 1431.487127 1196.286354  326.141932   30.607603 
#        utr3        utr5 
#   36.402326    8.601558 

# Write a BED file suitable for bedenrich
write.table(cbind(as.character(seqnames(final)), start(final), end(final), final$type),
    file=out.bed, sep='\t', col.names=F, row.names=F, quote=F)
