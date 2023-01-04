#!/usr/bin/env Rscript

# detect script being run by snakemake
# if so, make a mock commandArgs function
if ('snakemake' %in% ls()) {
    logfile <- snakemake@log[[1]]
    con <- file(logfile, 'w')
    sink(con, type='output')
    sink(con, type='message')

    commandArgs <- function(...) unlist(c(
        snakemake@input['neuron_histone'],
        snakemake@input['neuron_nott'],
        snakemake@input['neuron_scatacseq'],
        snakemake@input['neuron_scrnaseq'],
        snakemake@input['neuron_gtex'],
        snakemake@input['neuron_replichip'],
        snakemake@input['neuron_repliseq'],
        snakemake@input['oligo_histone'],
        snakemake@input['oligo_nott'],
        snakemake@input['oligo_scatacseq'],
        snakemake@input['oligo_scrnaseq'],
        snakemake@input['oligo_gtex'],
        snakemake@input['oligo_replichip'],
        snakemake@input['oligo_repliseq'],
        snakemake@output['pdf'],
        snakemake@output['svg']
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 16) {
    stop('usage: fig3_panel_e.R {neuron signals: histone, nott, scatacseq, scrnaseq, gtex, replichip, repliseq} {oligo signals: <same as neuron>} out.pdf out.svg')
}


neuron.fs <- args[1:7]
names(neuron.fs) <- c('histone', 'nott', 'scatacseq', 'scrnaseq', 'gtex', 'replichip', 'repliseq')
oligo.fs <- args[8:14]
names(oligo.fs) <- c('histone', 'nott', 'scatacseq', 'scrnaseq', 'gtex', 'replichip', 'repliseq')
out.pdf <- args[15]
out.svg <- args[16]

for (f in c(out.pdf, out.svg)) {
    if (file.exists(f))
        stop(paste('output file', f, 'already exists, please delete it first'))
}

suppressMessages(library(data.table))
suppressMessages(library(mutenrich))
suppressMessages(library(extrafont))
suppressMessages(library(svglite))
if (!("Arial" %in% fonts()))
    stop("Arial font not detected; did you load extrafonts and run font_import() with the appropriate path?")

library(data.table)

# BINSIZE = 1000   -- not doing 1000000 anymore
# ourcell = neuron, oligo
# IF <ourcell == neuron>
#     sigclass = SBS5, SBS12, SBS16, SBS89
# IF <ourcell == oligo>
#     sigclass = SBS1, SBS5, SBS32, SBS89
# atac_celltype = oligo, OPC, excitatory_neuron, inhibitory_neuron
# rna_celltype = OPCs, Oligodendrocytes, Excitatory-Neurons, Inhibitory-Neurons

# <BINSIZE> <sigclass>
# output format:
# roadmap_histone eid=E073 inactive(mark=H3K27me3 | mark=H3K9me3)
# roadmap_histone eid=E073 active(not above)
# nott celltype=<celltype> active(mark=H3K27ac | mark=H3K4me3)
# scatacseq libid=librarymerged celltype=<atac_celltype>
# nott mark=ATAC celltype=<celltype>
# scrnaseq celltype=<rna_celltype> letter(over donor,selection)
# gtex tissue=Brain_-_Frontal_Cortex__BA9_
# IF <BINSIZE == 1000000>
#    replichip neural_progenitor(encid=ENCFF469TYS | encid=ENCFF907YXU)
#        ENCFF203XOM -> ENCFF469TYS
#        ENCFF320GYL ->
#        ENCFF580VOP ->
#        ENCFF702XRT -> ENCFF907YXU

# set pdev=NULL to allow the caller to layout the panel
make.panel <- function(files, sig, ourcell=c('neuron', 'oligo'), pdev=x11, plot.1mb=TRUE) {
    ourcell <- match.arg(ourcell)
    col <- ifelse(ourcell == 'neuron', 1, 2)

    hst <- fread(files['histone'])[eid=='E073' & quantile %in% 1:10 & sigclass==sig][order(as.integer(quantile))][, c('mutsfrom', 'muttype') := list(ourcell, 'snv')]

    nott_celltype <- ourcell
    nott <- fread(files['nott'])[celltype == nott_celltype & quantile %in% 1:10 & sigclass==sig][order(as.integer(quantile))][, c('mutsfrom', 'muttype') := list(ourcell, 'snv')]

    atac_celltypes <- if (ourcell == 'oligo') c('OPC', 'oligo') else c('excitatory_neuron', 'inhibitory_neuron')
    atac <- fread(files['scatacseq'])[libid=='librarymerged' & celltype %in% atac_celltypes & quantile %in% 1:10 & sigclass==sig][order(as.integer(quantile))][, c('mutsfrom', 'muttype') := list(ourcell, 'snv')]

    rna_celltypes <- if (ourcell == 'oligo') c('OPCs', 'Oligodendrocytes') else c('Excitatory-Neurons', 'Inhibitory-Neurons')
    rna <- fread(files['scrnaseq'])[celltype %in% rna_celltypes & quantile %in% 1:10 & sigclass==sig][order(as.integer(quantile))][, c('mutsfrom', 'muttype') := list(ourcell, 'snv')]
    # Make a line for each donor,selection.celltype combo
    #rna$linegroup <- paste(rna$donor, rna$selection, rna$celltype)
    # Make an average line for each celltype
    rna <- rna[,.(enr=mean(enr)), by=.(BINSIZE, quantile, celltype)]
#print(rna)

    # GTEx using just BA9
    #gtex <- fread(paste0('gtex08_', ourcell, '_A_3quantiles.csv'))[tissue=='Brain_-_Frontal_Cortex__BA9_' & quantile %in% 1:10 & sigclass==sig][order(as.integer(quantile))]
    gtex <- fread(files['gtex'])[group=='Brain' & quantile %in% 1:10 & sigclass==sig][order(as.integer(quantile))][, c('mutsfrom', 'muttype') := list(ourcell, 'snv')]

    repli <- fread(files['replichip'])[encid %in% c('ENCFF469TYS', 'ENCFF907YXU') & quantile %in% 1:10 & sigclass==sig][order(as.integer(quantile))][, c('mutsfrom', 'muttype') := list(ourcell, 'snv')]
    # have to reverse quantile to make 1=early .. N=late
    repli$quantile <- (3:1)[as.integer(repli$quantile)]
    repli <- repli[order(as.integer(quantile))]

    repliseq <- fread(files['repliseq'])[celltype == 'SK-N-SH' & quantile %in% 1:10 & sigclass==sig][order(as.integer(quantile))][, c('mutsfrom', 'muttype') := list(ourcell, 'snv')]
    # have to reverse quantile to make 1=early .. N=late
    repliseq$quantile <- (3:1)[as.integer(repliseq$quantile)]

    if (!is.null(pdev)) {
        if (plot.1mb) {
            pdev(width=15,height=8, pointsize=5)
            layout(matrix(1:27, nrow=3, byrow=F), height=c(2,2,3))
        } else {
            pdev(width=15, height=5, pointsize=5)
            layout(matrix(1:18, nrow=2, byrow=F), height=c(2,3))
        }
    }
    par(mar=c(5,2,3,0.1))

    # inactive marks
    plotfn(hst[mark %in% c('H3K27me3', 'H3K9me3') & BINSIZE==1000],
            typecol='mark',
            labtype='number', xlab='Thirds', ylab='', main='Inactive marks',
            add.legend=!plot.1mb, col=col)
    if (plot.1mb) {
        plotfn(hst[mark %in% c('H3K27me3', 'H3K9me3') & BINSIZE==1000000],
                typecol='mark',
                labtype='number', xlab='Thirds', ylab='', main='Inactive marks',
                add.legend=T, col=col)
    }

    # active marks
    plotfn(hst[!(mark %in% c('H3K27me3', 'H3K9me3')) & BINSIZE==1000],
            typecol='mark',
            labtype='number', xlab='Thirds', ylab='', main='Active marks',
            add.legend=!plot.1mb, col=col)
    if (plot.1mb) {
        plotfn(hst[!(mark %in% c('H3K27me3', 'H3K9me3')) & BINSIZE==1000000],
                typecol='mark',
                labtype='number', xlab='Thirds', ylab='', main='Active marks',
                add.legend=T, col=col)
    }

    # Nott active marks
    plotfn(nott[mark %in% c('H3K27ac', 'H3K4me3') & BINSIZE==1000],
            typecol='mark',
            labtype='number', xlab='Thirds', ylab='', main='Nott active marks',
            add.legend=!plot.1mb, col=col)
    if (plot.1mb) {
        plotfn(nott[mark %in% c('H3K27ac', 'H3K4me3') & BINSIZE==1000000],
                typecol='mark',
                labtype='number', xlab='Thirds', ylab='', main='Nott active marks',
                add.legend=T, col=col)
    }

    # scATAC-seq
    plotfn(atac[BINSIZE==1000],
            typecol='celltype',
            labtype='number', xlab='Thirds', ylab='', main='scATAC-seq',
            add.legend=!plot.1mb, col=col)
    if (plot.1mb) {
        plotfn(atac[BINSIZE==1000000],
                typecol='celltype',
                labtype='number', xlab='Thirds', ylab='', main='scATAC-seq',
                add.legend=T, col=col)
    }

    # Nott ATAC
    plotfn(nott[mark == 'ATAC' & BINSIZE==1000],
            typecol='mark',
            labtype='number', xlab='Thirds', ylab='', main='Nott ATAC',
            add.legend=!plot.1mb, col=col)
    if (plot.1mb) {
        plotfn(nott[mark == 'ATAC' & BINSIZE==1000000],
                typecol='mark',
                labtype='number', xlab='Thirds', ylab='', main='Nott ATAC',
                add.legend=T, col=col)
    }

    # scRNAseq
    plotfn(rna[BINSIZE==1000],
            #typecol='linegroup',
            typecol='celltype',
            labtype='number', xlab='Thirds', ylab='', main='scRNA-seq',
            add.legend=!plot.1mb, col=col)
    if (plot.1mb) {
        plotfn(rna[BINSIZE==1000000],
                #typecol='linegroup',
                typecol='celltype',
                labtype='number', xlab='Thirds', ylab='', main='scRNA-seq',
                add.legend=T, col=col)
    }

    # GTEx
    plotfn(gtex[BINSIZE==1000],
            typecol='tissue',
            labtype='number', xlab='Thirds', ylab='', main='GTEx expression',
            add.legend=!plot.1mb, col=col)
    if (plot.1mb) {
        plotfn(gtex[BINSIZE==1000000],
                typecol='tissue',
                labtype='number', xlab='Thirds', ylab='', main='GTEx expression',
                add.legend=T, col=col)
    }

    # RepliChip
    # there is no binsize=1kb because the probes are not dense enough
    #plot(1, pch=NA, xlab='', ylab='', xaxt='n', yaxt='n', bty='n')
    # now always plot 1MB so we have some reference to compare RepliSeq to
        plotfn(repli[BINSIZE==1000000],
                typecol='encid',
                labtype='number', linetype='average',
                xlab='Thirds', ylab='', main='Replication timing (1 MB, RepliChIP)',
                add.legend=!plot.1mb, col=col)
        # above is linetype=average, meaning add.legend never happens
        plot(1, pch=NA, xlab='', ylab='', xaxt='n', yaxt='n', bty='n')
        legend('topleft', legend='Average', pch=20, cex=0.8)
    if (plot.1mb) {
        plotfn(repli[BINSIZE==1000000],
                typecol='encid',
                labtype='number', linetype='average',
                xlab='Thirds', ylab='', main='Replication timing',
                add.legend=T, col=col)
    }

    # RepliSeq has 1kb bins, probe density not an issue
    plotfn(repliseq[BINSIZE==1000],
            typecol='celltype',
            labtype='number', linetype='average',
            xlab='Thirds', ylab='', main='Replication timing (RepliSeq)',
            add.legend=!plot.1mb, col=col)
        # above is linetype=average, meaning add.legend never happens
        plot(1, pch=NA, xlab='', ylab='', xaxt='n', yaxt='n', bty='n')
        legend('topleft', legend='Average', pch=20, cex=0.8)
    if (plot.1mb) {
        plotfn(repli[BINSIZE==1000000],
                typecol='celltype',
                labtype='number', linetype='average',
                xlab='Thirds', ylab='', main='Replication timing (RepliSeq)',
                add.legend=T, col=col)
    }
}



plotfn <- function(n, typecol, linetype=c('separate', 'average'), labtype=c('point','number'), add.legend=TRUE, col=1, ...) {
    labtype <- match.arg(labtype)
    linetype <- match.arg(linetype)

    n <- n[quantile %in% 1:10]
    n$quantile <- as.integer(n$quantile)

    #ylim <- range(n$enr)*c(0.95,1.05)
    # Manually set uniform y-axis
    ylim <- c(0.4,2.25)
    xlim <- range(c(n$quantile))
    ## changes depending on track
    n$type <- n[[typecol]]  # only allows one column for now

    if (linetype == 'average') {
        n$quantile <- as.integer(n$quantile)

        plot(n[,.(mean(quantile), mean(enr)), by=quantile][,.(V1,V2)],
            type='b', lwd=2, col=col, pch=20, xlim=xlim, ylim=ylim, ...)
        abline(h=1, lty='dashed', col='grey')
    } else {
        types <- sort(unique(n$type))

        for (i in 1:length(types)) {
            t <- types[i]
            pf <- if (i == 1) plot else lines
            pch <- ifelse(labtype == 'point', 20, letters[i])
            pf(n[type == t, .(quantile, enr)], type='b', lwd=2, col=col, pch=pch,
                xlim=xlim, ylim=ylim, ...)
        }
        abline(h=1, lty='dashed', col='grey')

        if (add.legend) {
            # separate panel because it's so large
            plot(1, pch=NA, xlab='', ylab='', xaxt='n', yaxt='n', bty='n')
            legend('topleft', pch=letters[1:length(types)], legend=types, cex=0.8)
        }
    }
}


devs=list(pdf, svglite)
outs=c(out.pdf, out.svg)
for (i in 1:2) {
    devs[[i]](width=10, height=9, pointsize=5, file=outs[i])
    layout(rbind(
        matrix(1:18, nrow=2, byrow=F),
        18+matrix(1:18, nrow=2, byrow=F),
        36+matrix(1:18, nrow=2, byrow=F)), height=c(2,3,2,3,2,3))
    par(mar=c(5,2,3,0.1))
    # did more sigs for exploratory analysis, now just do panels for paper
    #for (sig in c('SBS1', 'SBS5', 'SBS32', 'SBS89')) {
    for (sig in c('SBS1', 'SBS32')) {
        print(c('oligo', sig))
        #make.panel(sig=sig, ourcell='oligo', pdev=function(...) pdf(file=paste0('oligo_', sig, '_v2.pdf'), ...))
        make.panel(files=oligo.fs, sig=sig, ourcell='oligo', pdev=NULL, plot.1mb=FALSE)
        #dev.off()
    }
    #for (sig in c('SBS5', 'SBS12', 'SBS16', 'SBS89')) {
    for (sig in 'SBS16') {
        print(c('neuron', sig))
        #make.panel(sig=sig, ourcell='neuron', pdev=function(...) pdf(file=paste0('neuron', sig, '_v2.pdf'), ...))
        make.panel(files=neuron.fs, sig=sig, ourcell='neuron', pdev=NULL, plot.1mb=FALSE)
        #dev.off()
    }

    dev.off()
}

if ('snakemake' %in% ls()) {
    sink()
}
