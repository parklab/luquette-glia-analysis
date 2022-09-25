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
        snakemake@params['extend_bp'],
        snakemake@output['csv'],
        snakemake@params['output_prefix']
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 5) {
    stop('usage: plot_gviz_region.R neuron_profile.csv oligo_profile.csv extend_basepairs out_regions.csv out_prefix_for_plots')
}


neuron.csv <- args[1]
oligo.csv <- args[2]
extend <- as.integer(args[3])
out.csv <- args[4]
out.prefix <- args[5]

for (f in c(out.csv)) {
    if (file.exists(f))
        stop(paste('output file', f, 'already exists, please delete it first'))
}

suppressMessages(library(Gviz))
suppressMessages(library(scan2))
suppressMessages(library(mutenrich))
suppressMessages(library(extrafont))
suppressMessages(library(svglite))
if (!("Arial" %in% fonts()))
    stop("Arial font not detected; did you load extrafonts and run font_import() with the appropriate path?")


get.annots <- function(chrom, from, to) {
    gtrack <- GenomeAxisTrack()
    itrack <- IdeogramTrack(genome='hg19', chromosome=chrom)
    refseq <- UcscTrack(genome='hg19', chromosome=chrom, from=from, to=to,
        table='ncbiRefSeq', track='NCBI RefSeq', trackType='GeneRegionTrack',
        rstarts='exonStarts', rends='exonEnds', gene='name2', symbol='name2',
        transcript='name', strand='strand', fill='#8282d2', name='RefSeq')
    
    list(itrack,gtrack,refseq)
}

ns <- fread(neuron.csv)
ns$chr <- paste0('chr', ns$chr)
os <- fread(oligo.csv)
os$chr <- paste0('chr', os$chr)
s <- ns[,.(chr,start,end,muts.smoothed.peryear)]
s$oligo.muts.smoothed.peryear <- os$muts.smoothed.peryear

dstrack <- DataTrack(GRanges(s), name='Mutations / Mb*year', col=1:2, groups=c('neuron','oligo'), type='l')

# Create full set of sites with Poisson enrichment p-value < 0.05 after BH adjustment
# TODO: make this configurable later. Allow selection of cutoff, but more importantly,
# whether to use one sided (enrichment only) or two sided (enrichment+depletion) tests.
ns[, CellType := 'Neuron']
os[, CellType := 'Oligo']
fwrite(rbind(ns,os)[pois.pval.enrich.adj < 0.05], file=out.csv)

all.sites <- reduce(GRanges(rbind(ns,os)[pois.pval.adj < 0.1]))

for (i in 1:length(all.sites)) {
    print(all.sites[i,])
    cat('fetching..\n')
    chrom <- as.character(seqnames(all.sites[i,]))
    from <- start(all.sites[i,])
    to <- end(all.sites[i,])
    cat('chrom', chrom, '\n')
    cat('from', from, '\n')
    cat('to', to, '\n')
    as <- get.annots(chrom=chrom, from=max(0, from-(2*extend)), to=to+extend*2)
    Sys.sleep(15)  # UCSC will block if too many requests are sent
    cat('plotting..\n')
    
    devs=list(pdf, svglite)
    outs=paste0(out.prefix, 'region_', i, c('.pdf', '.svg'))
    for (j in 1:2) {
        devs[[j]](width=4, height=6, pointsize=5, file=outs[j])
        #par(mar=c(4,4,3,1))
        ht <- HighlightTrack(trackList = list(dstrack, as[[3]]),
                        start = from, end=to, chromosome=chrom)
        plotTracks(c(as[1:2],ht),
            from=max(0, from-extend), to=to+extend,
            collapseTranscripts='meta',
            transcriptAnnotation='symbol',
            legend=FALSE)
        dev.off()
    }
}

if ('snakemake' %in% ls()) {
    sink()
}
