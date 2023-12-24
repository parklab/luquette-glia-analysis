#!/usr/bin/env Rscript

# detect script being run by snakemake
# if so, make a mock commandArgs function
if ('snakemake' %in% ls()) {
    logfile <- snakemake@log[[1]]
    con <- file(logfile, 'w')
    sink(con, type='output')
    sink(con, type='message')

    commandArgs <- function(...) unlist(c(
        snakemake@input['neuron_scrnaseq'],
        snakemake@input['neuron_scatacseq'],
        snakemake@input['neuron_repliseq'],
        snakemake@input['neuron_histone'],
        snakemake@input['oligo_scrnaseq'],
        snakemake@input['oligo_scatacseq'],
        snakemake@input['oligo_repliseq'],
        snakemake@input['oligo_histone'],
        snakemake@output['pdf'],
        snakemake@output['svg'],
        snakemake@output['csv']
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 11) {
    stop('usage: fig5_allpanels.R {neuron signals: scrnaseq, scatacseq, repliseq, histone} {oligo signals: <same as neuron>} out.pdf out.svg out.csv')
}


neuron.fs <- args[1:4]
names(neuron.fs) <- c('scrnaseq', 'scatacseq', 'repliseq', 'histone')
oligo.fs <- args[5:8]
names(oligo.fs) <- c('scrnaseq', 'scatacseq', 'repliseq', 'histone')
out.pdf <- args[9]
out.svg <- args[10]
out.csv <- args[11]

for (f in c(out.pdf, out.svg, out.csv)) {
    if (file.exists(f))
        stop(paste('output file', f, 'already exists, please delete it first'))
}

suppressMessages(library(data.table))
suppressMessages(library(mutenrich))
suppressMessages(library(svglite))

log.y.axis <- 'y'
#log.y.axis <- ''

# this doesn't actually work, just introduces some instrumentation. always uses
# elements 1 and 2 regardless of what groups are being plotted.
lty.map <- c(
`pta_neuron`='solid',
`pta_oligo`='dashed'
)

col.map <- c(
`Astrocytes`='#F4C013',
`Endothelial`='#65B648',
`Excitatory-Neurons`='#000000',
`Inhibitory-Neurons`='#52C0CC',
`Microglia`='#A8388C',
`OPCs`='#929393',
`Oligodendrocytes`='#E91D21',
`astrocyte`='#F4C013',
`endothelial`='#65B648',
`excitatory_neuron`='#000000',
`inhibitory_neuron`='#52C0CC',
`microglia`='#A8388C',
`OPC`='#929393',
`oligo`='#E91D21',
`BG02ES`="black",
`BJ`="black",
`GM06990`="black",
`GM12801`="black",
`GM12812`="black",
`GM12813`="black",
`GM12878`="black",
`HUVEC`="black",
`HeLa-S3`="black",
`HepG2`="black",
`IMR90`="black",
`K562`="black",
`MCF-7`="black",
`NHEK`="black",
`SK-N-SH`="black",
`H3K27ac`='#66c2a5',
`H3K36me3`='#fc8d62',
`H3K4me1`='#8da0cb',
`H3K4me3`='#e78ac3',
`H3K9ac`='#a6d854',
`H3K27me3`='#000000',
`H3K9me3`='#e5c494'
)


# Read in tables from various genomic covariates. These can have different
# metadata columns, so normalize them to remove metadata columns except:
#   datasource - one of: scrnaseq, scatacseq, repliseq, inactive_histone_mark, active_histone_mark
#   lineclass - for enrichment plots with multiple lines, this column defines each line
# 'ourcell' adds a column to the table recording whether these enrichment
#    analyses come from neurons or oligos.
consolidate.tables <- function(files, ourcell=c('neuron', 'oligo')) {
    ourcell <- match.arg(ourcell)

    wrap.fread <- function(...) {
        fread(...)[quantile %in% 1:100][order(as.integer(quantile))][, c('mutsfrom', 'muttype') := list(ourcell, 'snv')][BINSIZE == 1000]
    }

    hst <- wrap.fread(files['histone'])[eid == 'E073']
    hst$lineclass <- hst$mark
    hst$datasource <- ifelse(hst$mark %in% c('H3K27me3', 'H3K9me3'), 'inactive_histone_mark', 'active_histone_mark')
    hst$dataclass <- NULL
    hst$signal_type <- NULL
    hst$eid <- NULL
    hst$mark <- NULL

    atac <- wrap.fread(files['scatacseq'])[libid == 'librarymerged']
    map <- c(OPC='OPCs', astrocyte='Astrocytes', endothelial='Endothelial',
        excitatory_neuron='Excitatory-Neurons', inhibitory_neuron='Inhibitory-Neurons', 
        microglia="Microglia", oligo="Oligodendrocytes")
    atac$celltype <- map[atac$celltype]
    atac$lineclass <- atac$celltype
    atac$libid <- NULL
    atac$sample <- NULL
    atac$celltype <- NULL

    # scRNAseq
    rna <- wrap.fread(files['scrnaseq'])[donor == 'combined' & selection == 'combined']
    rna$lineclass <- rna$celltype
    rna$dataclass <- NULL
    rna$donor <- NULL
    rna$selection <- NULL
    rna$celltype <- NULL


    # RepliSeq from ENCODE. 
    repliseq <- wrap.fread(files['repliseq'])#[celltype == 'SK-N-SH']  # shouldn't use SK-N-SH because fig4 doesn't
    # have to reverse quantile to make 1=early .. N=late
    repliseq$quantile <- (3:1)[as.integer(repliseq$quantile)]
    repliseq$lineclass <- paste(repliseq$celltype, repliseq$geoid)
    repliseq$celltype <- NULL
    repliseq$geoid <- NULL

    rbind(rna, atac, repliseq, hst)
}


make.panels <- function(tab, sig, colors=c(neuron=1, oligo=2), add.legend=FALSE, show.title=FALSE, show.xaxis=FALSE) {
print(sig)
print(tab[sigclass == sig & datasource != 'repliseq']$enr)
    # repliseq has one major outlier BEFORE averaging. it never defines the windows
    # in these data, so as a hack just don't use it in getting ylims.
    #ylims <- range(tab[sigclass == sig & datasource != 'repliseq']$enr)*c(0.95,1.05)
    ylims <- range(pretty(tab[sigclass == sig & datasource != 'repliseq']$enr)) #*c(0.9,1.3)
    if (ylims[1] == 0 & log.y.axis == 'y')
        ylims[1] <- 0.3
print(ylims)

    # scRNAseq
    plotfn(tab[datasource == 'scrnaseq' & sigclass == sig],
            labtype='number', main='scRNA-seq', show.title=show.title,
            add.legend=add.legend, ylim=ylims, ylab='Obs / Exp',
            yaxis=TRUE, xaxis=show.xaxis)

    # scATAC-seq
    plotfn(tab[datasource == 'scatacseq' & sigclass == sig],
            labtype='number', main='scATAC-seq',
            add.legend=add.legend, ylim=ylims, show.title=show.title,
            ylab='', yaxis=FALSE, xaxis=show.xaxis)

    # RepliSeq
    plotfn(tab[datasource == 'repliseq' & sigclass == sig],
            labtype='number', linetype='average',
            main='Replication timing', show.title=show.title,
            add.legend=add.legend, ylim=ylims,
            ylab='', yaxis=FALSE, xaxis=show.xaxis)

    # active marks
    plotfn(tab[datasource == 'active_histone_mark' & sigclass == sig],
            labtype='number', main='Active histone marks',
            add.legend=add.legend, ylim=ylims, show.title=show.title,
            ylab='', yaxis=FALSE, xaxis=show.xaxis)

    # inactive marks
    plotfn(tab[datasource == 'inactive_histone_mark' & sigclass == sig],
            labtype='number', main='Inactive histone marks',
            add.legend=add.legend, ylim=ylims, show.title=show.title,
            ylab='', yaxis=FALSE, xaxis=show.xaxis)
}



# always using lty.map[1]/[2]  doesn't really allow for flexible plotting, it just
# introduces some of the instrumentation required for that eventually.
plotfn <- function(n, ylim, xlab='', ylab='', main='', typecol='lineclass', linetype=c('separate', 'average'), labtype=c('point','number'), add.legend=TRUE, col=1, yaxis=FALSE, xaxis=FALSE, neuron.col=1, oligo.col=2, show.title=FALSE, yaxt='s', xaxt='s', add=FALSE, lwd=2/3, group1.lty=lty.map[1], group2.lty=lty.map[2], ...) {
    labtype <- match.arg(labtype)
    linetype <- match.arg(linetype)

    xlim <- range(c(n$quantile))
    ## changes depending on track
    n$type <- n[[typecol]]  # only allows one column for now

    if (show.title == FALSE)
        main <- ''

    if (xaxis == FALSE) {
        xaxt <- 'n'
        xlab <- ''
    }
    if (yaxis == FALSE) {
        yaxt <- 'n'
        ylab <- ''
    }

    make.y.axis <- FALSE
    if (yaxis == TRUE & log.y.axis == 'y') {
        make.y.axis <- TRUE
        yaxt <- 'n'
    }
        

    if (show.title) {
        # separate panel because it's so large
        par(mar=c(0,0,0,0))
        plot(1, pch=NA, xlab='', ylab='', xaxt='n', yaxt='n', bty='n')
        legend('center', legend=main, bty='n')
    }
    main <- ''

    par(mar=c(2, 4, 0.1, 0))
    if (linetype == 'average') {
        n$quantile <- as.integer(n$quantile)

        plotfn <- plot
        if (add == TRUE)
            plotfn <- lines

        # averaging over many dataclasses (=type column), arbitrarily use the color of the first
        t <- strsplit(n$type[1], ' ')[[1]][1]
        group1.col <- col.map[t]
        group2.col <- col.map[t]
cat('average line: selecting color for type t=', t, 'color=', group1.col, '\n')

        plotfn(n[mutsfrom == 'neuron',.(mean(quantile), mean(enr)), by=quantile][,.(V1,V2)],
            type='l', lwd=lwd, col=group1.col, pch=20, xlim=xlim, ylim=ylim, bty='n', xlab=xlab, ylab=ylab, main=main, xaxt=xaxt, yaxt=yaxt, lty=group1.lty, log=log.y.axis, ...)
        if (nrow(n[mutsfrom == 'oligo']) > 0) {
            lines(n[mutsfrom == 'oligo',.(mean(quantile), mean(enr)), by=quantile][,.(V1,V2)],
                type='l', lwd=lwd, col=group2.col, pch=20, xlim=xlim, ylim=ylim, bty='n', xlab=xlab, ylab=ylab, main=main, xaxt=xaxt, yaxt=yaxt, lty=group2.lty, log=log.y.axis, ...)
        } else {
            cat('skipping mutsfrom=oligo, nrow=0\n')
        }
        abline(h=1, lty='solid', col='black', lwd=1/2)
        if (make.y.axis) {
            axis(side=2, at=pretty(ylim, n=7))
        }
        # make a blank legend just to honor the expectation of add.legend=TRUE
        if (add.legend) {
            # separate panel because it's so large
            plot(1, pch=NA, xlab='', ylab='', xaxt='n', yaxt='n', bty='n')
        }
    } else {
        types <- sort(unique(n$type))

        for (i in 1:length(types)) {
            t <- types[i]
            plot.type <- 'b'
            group1.col <- col.map[t]
            group2.col <- col.map[t]
            plot.type <- 'l'   # never use the letters anymore
cat('plotting trace for type t=', t, 'color=', group1.col, '\n')
            pf <- if (i == 1 & !add) plot else lines
            pch <- ifelse(labtype == 'point', 20, letters[i])
cat('neuron -----------------------------------------------\n')
print(n[mutsfrom == 'neuron' & type == t, .(quantile, enr)])
cat('oligo  -----------------------------------------------\n')
print(n[mutsfrom == 'oligo' & type == t, .(quantile, enr)])
            pf(n[mutsfrom == 'neuron' & type == t, .(quantile, enr)], type=plot.type, lwd=lwd, col=group1.col, pch=pch,
                xlim=xlim, ylim=ylim, bty='n', xlab=xlab, ylab=ylab, main=main, xaxt=xaxt, yaxt=yaxt, lty=group1.lty, log=log.y.axis, ...)
            if (nrow(n[mutsfrom == 'oligo' & type == t]) > 0) {
                lines(n[mutsfrom == 'oligo' & type == t, .(quantile, enr)], type=plot.type, lwd=lwd, col=group2.col, pch=pch,
                    xlim=xlim, ylim=ylim, bty='n', xlab=xlab, ylab=ylab, main=main, xaxt=xaxt, yaxt=yaxt, lty=group2.lty, log=log.y.axis, ...)
            } else {
                cat('skipping mutsfrom=oligo type=', t, ', nrow=0\n')
            }
        }
        abline(h=1, lty='solid', col='black', lwd=1/2)
        if (make.y.axis) {
            axis(side=2, at=pretty(ylim, n=7))
        }

        if (add.legend) {
            # separate panel because it's so large
            par(mar=c(0,0,0,0))
            plot(1, pch=NA, xlab='', ylab='', xaxt='n', yaxt='n', bty='n')
            legend('center', fill=col.map[types], legend=types, border=col.map[types], bty='n')
            #legend('topleft', pch=letters[1:length(types)], legend=types, bty='n')
        }
    }
}

sigs.to.plot <- c('SBS1', 'SBS16', 'SBS5')

n <- consolidate.tables(files=neuron.fs, ourcell='neuron')
o <- consolidate.tables(files=oligo.fs, ourcell='oligo')
all.data <- rbind(n, o)[lineclass != "Neurons"]  # there was one scRNAseq library that did not differentiate excitatory and inhibitory neurons

#stop('done')

fwrite(all.data, file=out.csv)

devs=list(pdf, svglite)
outs=c(out.pdf, out.svg)
for (i in 1:2) {
    devs[[i]](width=6.5, height=(6.5/5)*4+0.1, pointsize=5, file=outs[i])

    layout(rbind(matrix(1:10, nrow=2), 10+1:5, 15+matrix(1:10,nrow=2)), height=c(1,5,5,5,5))
    #make.panels(tab=all.data[mutsfrom=='oligo'], sig=sigs.to.plot[1], show.title=TRUE, show.xaxis=FALSE)
    make.panels(tab=all.data, sig=sigs.to.plot[1], show.title=TRUE, show.xaxis=FALSE)
    make.panels(tab=all.data[mutsfrom=='neuron'], sig=sigs.to.plot[2], show.title=FALSE, show.xaxis=FALSE)
    #make.panels(tab=all.data, sig=sigs.to.plot[2], show.title=FALSE, show.xaxis=FALSE)
    make.panels(tab=all.data, sig=sigs.to.plot[3], show.title=FALSE, show.xaxis=TRUE, add.legend=TRUE)
    dev.off()
}

warnings()

if ('snakemake' %in% ls()) {
    sink()
}
