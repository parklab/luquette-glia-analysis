library(Gviz)
library(scan2)
library(mutenrich)

extend=1e6

get.annots <- function(chrom, from, to) {
    gtrack <- GenomeAxisTrack()
    itrack <- IdeogramTrack(genome='hg19', chromosome=chrom)
    refseq <- UcscTrack(genome='hg19', chromosome=chrom, # from=from, to=to,
        table='ncbiRefSeq', track='NCBI RefSeq', trackType='GeneRegionTrack',
        rstarts='exonStarts', rends='exonEnds', gene='name2', symbol='name2',
        transcript='name', strand='strand', fill='#8282d2', name='RefSeq')
    
    list(itrack,gtrack,refseq)
}

ns <- fread('neuron.profile.csv')
ns$chr <- paste0('chr', ns$chr)
os <- fread('oligo.profile.csv')
os$chr <- paste0('chr', os$chr)
s <- ns[,.(chr,start,end,muts.smoothed.peryear)]
s$oligo.muts.smoothed.peryear <- os$muts.smoothed.peryear

ni <- fread('neuron_indel_passAB.profile.csv')
ni$chr <- paste0('chr', ni$chr)
oi <- fread('oligo_indel_passAB.profile.csv')
oi$chr <- paste0('chr', oi$chr)
i <- ni[,.(chr,start,end,muts.smoothed.peryear)]
i$oligo.muts.smoothed.peryear <- oi$muts.smoothed.peryear


dstrack <- DataTrack(GRanges(s), name='SNVs per Mb*year', col=1:2, groups=c('neuron','oligo'), type='l')
ditrack <- DataTrack(GRanges(i), name='Indels per Mb*year', col=1:2, groups=c('neuron','oligo'), type='l')

# Create full set of sites with Poisson p-value < 0.1 after BH adjustment
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
    as <- get.annots(chrom=chrom)
    cat('plotting..\n')

    ht <- HighlightTrack(trackList = list(dstrack, ditrack, as[[3]]),
                     start = from, end=to, chromosome=chrom)
    pdf(file=paste0('region_', i, '.pdf'))
    plotTracks(c(as[1:2],ht),
        from=from-extend, to=to+extend,
        collapseTranscripts='meta',
        transcriptAnnotation='symbol',
        legend=FALSE)
    dev.off()
}
