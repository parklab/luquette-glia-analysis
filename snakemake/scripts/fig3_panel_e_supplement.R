#!/usr/bin/env Rscript

# detect script being run by snakemake
# if so, make a mock commandArgs function
if ('snakemake' %in% ls()) {
    logfile <- snakemake@log[[1]]
    con <- file(logfile, 'w')
    sink(con, type='output')
    sink(con, type='message')

    commandArgs <- function(...) unlist(c(
        snakemake@params['celltype'],
        snakemake@params['signame'],
        snakemake@input['histone'],
        snakemake@input['nott'],
        snakemake@input['scatacseq'],
        snakemake@input['scrnaseq'],
        snakemake@input['gtex'],
        snakemake@input['replichip'],
        snakemake@input['repliseq'],
        snakemake@output['pdf'],
        snakemake@output['svg']
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 11) {
    stop('usage: fig3_panel_e_supplement.R {neuron|oligo} sig_name histone.csv nott.csv scatacseq.csv scrnaseq.csv gtex.csv replichip.csv repliseq.csv out.pdf out.svg')
}


celltype <- args[1]
signame <- args[2]
in.fs <- args[3:9]
names(in.fs) <- c('histone', 'nott', 'scatacseq', 'scrnaseq', 'gtex', 'replichip', 'repliseq')
out.pdf <- args[10]
out.svg <- args[11]

for (f in c(out.pdf, out.svg)) {
    if (file.exists(f))
        stop(paste('output file', f, 'already exists, please delete it first'))
}

suppressMessages(library(data.table))
suppressMessages(library(mutenrich))
suppressMessages(library(svglite))

library(data.table)

# atac_celltype = oligo, OPC, excitatory_neuron, inhibitory_neuron
# rna_celltype = OPCs, Oligodendrocytes, Excitatory-Neurons, Inhibitory-Neurons


# Tracks we use:
# roadmap_histone eid=E073 inactive(mark=H3K27me3 | mark=H3K9me3)
# roadmap_histone eid=E073 active(not above)
# nott celltype=<celltype> active(mark=H3K27ac | mark=H3K4me3)
# scatacseq libid=librarymerged celltype=<atac_celltype>
# nott mark=ATAC celltype=<celltype>
# scrnaseq celltype=<rna_celltype> letter(over donor,selection)
# gtex tissue=Brain_-_Frontal_Cortex__BA9_

# set pdev=NULL to allow the caller to layout the panel
#
# all.tracks - plot all tracks, not the suggested tracks above
make.panel <- function(files, sig, ourcell=c('neuron', 'oligo'), pdev=x11, plot.1mb=TRUE, all.tracks=TRUE) {
    ourcell <- match.arg(ourcell)
    col <- ifelse(ourcell == 'neuron', 1, 2)

    # This code is stolen from the sigenrich analysis in which the tables of
    # enrichment levels have a "sigclass" column corresponding to each signature in the
    # reduced COSMIC set. This script plots the output of a regular {,q}bedenrich that is
    # run on {celltype}_sbs{sig} mutation sets, which are just a subset of the complete
    # mutation calls that approximate a particular signature (e.g., C>T at CpGs for SBS1).
    # These tables do not have the requisite sigclass column, so add it here.
    wrap.fread <- function(...) {
        fread(...)[, sigclass := sig][quantile %in% 1:50][order(as.integer(quantile))][, c('mutsfrom', 'muttype') := list(ourcell, 'snv')]
    }

    hst <- wrap.fread(files['histone'])
    if (!all.tracks)
        hst <- hst[eid == 'E073']

    nott_celltype <- ourcell
    nott <- wrap.fread(files['nott'])
    if (!all.tracks)
        nott <- nott[celltype == nott_celltype]

    atac_celltypes <- if (ourcell == 'oligo') c('OPC', 'oligo') else c('excitatory_neuron', 'inhibitory_neuron')
    atac <- wrap.fread(files['scatacseq'])[libid=='librarymerged']
    if (!all.tracks)
        atac <- atac[celltype %in% atac_celltypes]

    rna_celltypes <- if (ourcell == 'oligo') c('OPCs', 'Oligodendrocytes') else c('Excitatory-Neurons', 'Inhibitory-Neurons')
    rna <- wrap.fread(files['scrnaseq'])
    if (!all.tracks)
        rna <- rna[celltype %in% rna_celltypes]
    # Make a line for each donor,selection.celltype combo
    #rna$linegroup <- paste(rna$donor, rna$selection, rna$celltype)
    # Make an average line for each celltype
    rna <- rna[,.(enr=mean(enr)), by=.(BINSIZE, quantile, celltype)]

    # GTEx using just BA9
    gtex <- wrap.fread(files['gtex'])[group=='Brain']
    if (!all.tracks)
        gtex <- gtex[tissue=='Brain_-_Frontal_Cortex__BA9_']

    repli <- wrap.fread(files['replichip'])
    # have to reverse quantile to make 1=early .. N=late
    repli$quantile <- (max(repli$quantile):1)[as.integer(repli$quantile)]
    repli <- repli[order(as.integer(quantile))]
    if (!all.tracks)
        repli <- repli[encid %in% c('ENCFF469TYS', 'ENCFF907YXU')]

    repliseq <- wrap.fread(files['repliseq'])
    # have to reverse quantile to make 1=early .. N=late
    repliseq$quantile <- (max(repliseq$quantile):1)[as.integer(repliseq$quantile)]
    if (!all.tracks)
        repliseq <- repliseq[celltype == 'SK-N-SH']

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

    xlab <- function(data) {
        nq <- length(unique(data$quantile))
print(unique(data$quantile))
print(nq)
        if (nq == 3) "Thirds"
        else if (nq == 4) "Quartiles"
        else if (nq == 5) "Pentiles"
        else if (nq == 10) "Deciles"
        else paste0(as.character(nq), 'iles')
    }

    # inactive marks
    plotfn(hst[mark %in% c('H3K27me3', 'H3K9me3') & BINSIZE==1000],
            typecol='mark',
            labtype='number', xlab=xlab(hst), ylab='', main='Inactive marks',
            add.legend=!plot.1mb, col=col)
    if (plot.1mb) {
        plotfn(hst[mark %in% c('H3K27me3', 'H3K9me3') & BINSIZE==1000000],
                typecol='mark',
                labtype='number', xlab=xlab(hst), ylab='', main='Inactive marks',
                add.legend=T, col=col)
    }

    # active marks
    plotfn(hst[!(mark %in% c('H3K27me3', 'H3K9me3')) & BINSIZE==1000],
            typecol='mark',
            labtype='number', xlab=xlab(hst), ylab='', main='Active marks',
            add.legend=!plot.1mb, col=col)
    if (plot.1mb) {
        plotfn(hst[!(mark %in% c('H3K27me3', 'H3K9me3')) & BINSIZE==1000000],
                typecol='mark',
                labtype='number', xlab=xlab(hst), ylab='', main='Active marks',
                add.legend=T, col=col)
    }

    # Nott active marks
    plotfn(nott[mark %in% c('H3K27ac', 'H3K4me3') & BINSIZE==1000],
            typecol='mark',
            labtype='number', xlab=xlab(nott), ylab='', main='Nott active marks',
            add.legend=!plot.1mb, col=col)
    if (plot.1mb) {
        plotfn(nott[mark %in% c('H3K27ac', 'H3K4me3') & BINSIZE==1000000],
                typecol='mark',
                labtype='number', xlab=xlab(nott), ylab='', main='Nott active marks',
                add.legend=T, col=col)
    }

    # scATAC-seq
    plotfn(atac[BINSIZE==1000],
            typecol='celltype',
            labtype='number', xlab=xlab(atac), ylab='', main='scATAC-seq',
            add.legend=!plot.1mb, col=col)
    if (plot.1mb) {
        plotfn(atac[BINSIZE==1000000],
                typecol='celltype',
                labtype='number', xlab=xlab(atac), ylab='', main='scATAC-seq',
                add.legend=T, col=col)
    }

    # Nott ATAC
    plotfn(nott[mark == 'ATAC' & BINSIZE==1000],
            typecol='mark',
            labtype='number', xlab=xlab(nott), ylab='', main='Nott ATAC',
            add.legend=!plot.1mb, col=col)
    if (plot.1mb) {
        plotfn(nott[mark == 'ATAC' & BINSIZE==1000000],
                typecol='mark',
                labtype='number', xlab=xlab(nott), ylab='', main='Nott ATAC',
                add.legend=T, col=col)
    }

    # scRNAseq
    plotfn(rna[BINSIZE==1000],
            #typecol='linegroup',
            typecol='celltype',
            labtype='number', xlab=xlab(rna), ylab='', main='scRNA-seq',
            add.legend=!plot.1mb, col=col)
    if (plot.1mb) {
        plotfn(rna[BINSIZE==1000000],
                #typecol='linegroup',
                typecol='celltype',
                labtype='number', xlab=xlab(rna), ylab='', main='scRNA-seq',
                add.legend=T, col=col)
    }

    # GTEx
    plotfn(gtex[BINSIZE==1000],
            typecol='tissue',
            labtype='number', xlab=xlab(gtex), ylab='', main='GTEx expression',
            add.legend=!plot.1mb, col=col)
    if (plot.1mb) {
        plotfn(gtex[BINSIZE==1000000],
                typecol='tissue',
                labtype='number', xlab=xlab(gtex), ylab='', main='GTEx expression',
                add.legend=T, col=col)
    }

    # RepliChip
    # there is no binsize=1kb because the probes are not dense enough
    #plot(1, pch=NA, xlab='', ylab='', xaxt='n', yaxt='n', bty='n')
    # now always plot 1MB so we have some reference to compare RepliSeq to
        plotfn(repli[BINSIZE==1000000],
                typecol='encid',
                labtype='number', linetype='average',
                xlab=xlab(repli), ylab='', main='Replication timing (1 MB, RepliChIP)',
                add.legend=!plot.1mb, col=col)
        # above is linetype=average, meaning add.legend never happens
        plot(1, pch=NA, xlab='', ylab='', xaxt='n', yaxt='n', bty='n')
        legend('topleft', legend='Average', pch=20, cex=0.8)
    if (plot.1mb) {
        plotfn(repli[BINSIZE==1000000],
                typecol='encid',
                labtype='number', linetype='average',
                xlab=xlab(repli), ylab='', main='Replication timing',
                add.legend=T, col=col)
    }

    # RepliSeq has 1kb bins, probe density not an issue
    plotfn(repliseq[BINSIZE==1000],
            typecol='celltype',
            labtype='number', linetype='average',
            xlab=xlab(repliseq), ylab='', main='Replication timing (RepliSeq)',
            add.legend=!plot.1mb, col=col)
        # above is linetype=average, meaning add.legend never happens
        plot(1, pch=NA, xlab='', ylab='', xaxt='n', yaxt='n', bty='n')
        legend('topleft', legend='Average', pch=20, cex=0.8)
    if (plot.1mb) {
        plotfn(repli[BINSIZE==1000000],
                typecol='celltype',
                labtype='number', linetype='average',
                xlab=xlab(repli), ylab='', main='Replication timing (RepliSeq)',
                add.legend=T, col=col)
    }
}



plotfn <- function(n, typecol, linetype=c('separate', 'average'), labtype=c('point','number'), add.legend=TRUE, col=1, ...) {
    labtype <- match.arg(labtype)
    linetype <- match.arg(linetype)

    n <- n[quantile %in% 1:10]
    n$quantile <- as.integer(n$quantile)

    ylim <- range(n$enr)*c(0.95,1.05)
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
    devs[[i]](width=10, height=3, pointsize=5, file=outs[i])
    #layout(rbind(
        #matrix(1:9, nrow=3, byrow=F),
        #9+matrix(1:9, nrow=3, byrow=F)),
        #height=c(2,3))
    layout(matrix(1:18, byrow=F, nrow=2), height=c(2,3))
    par(mar=c(5,2,3,0.1))

    make.panel(files=in.fs, sig=signame, ourcell=celltype, pdev=NULL, plot.1mb=FALSE, all.tracks=FALSE)

    dev.off()
}

if ('snakemake' %in% ls()) {
    sink()
}
