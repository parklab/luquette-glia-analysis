#!/usr/bin/env Rscript

# detect script being run by snakemake
# if so, make a mock commandArgs function
if ('snakemake' %in% ls()) {
    logfile <- snakemake@log[[1]]
    con <- file(logfile, 'w')
    sink(con, type='output')
    sink(con, type='message')

    commandArgs <- function(...) unlist(c(
        snakemake@input['neuron'],
        snakemake@input['oligo'],
        snakemake@output['csv'],
        snakemake@output['pdf'],
        snakemake@output['svg']
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 7) {
    stop('usage: plot_gviz_region.R neuron_snv_profile.csv neuron_indel_profile.csv oligo_snv_profile.csv oligo_indel_profile.csv extend_basepairs out_regions.csv out_prefix_for_plots')
}


neuron_snv.csv <- args[1]
neuron_indel.csv <- args[2]
oligo_snv.csv <- args[3]
oligo_indel.csv <- args[4]
extend <- as.integer(args[5])
out.csv <- args[6]
out.prefix <- args[7]

for (f in c(out.csv)) {
    if (file.exists(f))
        stop(paste('output file', f, 'already exists, please delete it first'))
}

suppressMessages(library(Gviz))
suppressMessages(library(scan2))
suppressMessages(library(mutenrich))
suppressMessages(library(svglite))


get.annots <- function(chrom, from, to) {
    gtrack <- GenomeAxisTrack()
    itrack <- IdeogramTrack(genome='hg19', chromosome=chrom)
    refseq <- UcscTrack(genome='hg19', chromosome=chrom, from=from, to=to,
        table='ncbiRefSeq', track='NCBI RefSeq', trackType='GeneRegionTrack',
        rstarts='exonStarts', rends='exonEnds', gene='name2', symbol='name2',
        transcript='name', strand='strand', fill='#8282d2', name='RefSeq')
    
    list(itrack,gtrack,refseq)
}

ns <- fread(neuron_snv.csv)
ns$chr <- paste0('chr', ns$chr)
os <- fread(oligo_snv.csv)
os$chr <- paste0('chr', os$chr)
s <- ns[,.(chr,start,end,muts.smoothed.peryear)]
s$oligo.muts.smoothed.peryear <- os$muts.smoothed.peryear

ni <- fread(neuron_indel.csv)
ni$chr <- paste0('chr', ni$chr)
oi <- fread(oligo_indel.csv)
oi$chr <- paste0('chr', oi$chr)
i <- ni[,.(chr,start,end,muts.smoothed.peryear)]
i$oligo.muts.smoothed.peryear <- oi$muts.smoothed.peryear


dstrack <- DataTrack(GRanges(s), name='SNVs per Mb*year', col=1:2, groups=c('neuron','oligo'), type='l')
ditrack <- DataTrack(GRanges(i), name='Indels per Mb*year', col=1:2, groups=c('neuron','oligo'), type='l')

# Create full set of sites with Poisson p-value < 0.1 after BH adjustment
fwrite(rbind(ns,os,ni,oi)[pois.pval.adj < 0.1], file=out.csv)

all.sites <- reduce(GRanges(rbind(ns,os,ni,oi)[pois.pval.adj < 0.1]))

for (i in 1:length(all.sites)) {
    print(all.sites[i,])
    cat('fetching..\n')
    chrom <- as.character(seqnames(all.sites[i,]))
    from <- start(all.sites[i,])
    to <- end(all.sites[i,])
    cat('chrom', chrom, '\n')
    cat('from', from, '\n')
    cat('to', to, '\n')
    as <- get.annots(chrom=chrom, from=from-(2*extend), to=to+extend*2)
    cat('plotting..\n')
    
    devs=list(pdf, svglite)
    outs=paste0(out.prefix, 'region_', i, c('.pdf', '.svg'))
    for (j in 1:2) {
        devs[[j]](width=4, height=6, pointsize=5, file=outs[j])
        #par(mar=c(4,4,3,1))
        ht <- HighlightTrack(trackList = list(dstrack, ditrack, as[[3]]),
                        start = from, end=to, chromosome=chrom)
        plotTracks(c(as[1:2],ht),
            from=from-extend, to=to+extend,
            collapseTranscripts='meta',
            transcriptAnnotation='symbol',
            legend=FALSE)
        dev.off()
    }
}

if ('snakemake' %in% ls()) {
    sink()
}
