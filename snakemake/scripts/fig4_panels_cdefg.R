#!/usr/bin/env Rscript

# detect script being run by snakemake
# if so, make a mock commandArgs function
if ('snakemake' %in% ls()) {
    logfile <- snakemake@log[[1]]
    con <- file(logfile, 'w')
    sink(con, type='output')
    sink(con, type='message')

    commandArgs <- function(...) unlist(c(
        snakemake@input['group1'],
        snakemake@input['group2'],
        snakemake@params['group1_color'],
        snakemake@params['group2_color'],
        snakemake@params['log_y_axis'],
        snakemake@output['pdf'],
        snakemake@output['svg'],
        snakemake@output['csv']
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 8) {
    cat('use group1_color=col_map to instead use a different color for each line and dashed/solid lines for OLs and neurons.\n')
    stop('usage: fig4_panels_cdefg.R group1_consolidated_table.csv group2_consolidated_table.csv group1_color group2_color {log_y_axis=true|false} out.pdf out.svg out.csv')
}


group1.csv <- args[1]
group2.csv <- args[2]
group1.col <- unname(args[3])
group2.col <- unname(args[4])
log.y.axis <- ifelse(as.logical(args[5]), 'y', '')
out.pdf <- args[6]
out.svg <- args[7]
out.csv <- args[8]

for (f in c(out.pdf, out.svg, out.csv)) {
    if (file.exists(f))
        stop(paste('output file', f, 'already exists, please delete it first'))
}

suppressMessages(library(data.table))
suppressMessages(library(mutenrich))
suppressMessages(library(svglite))


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
`BG02ES`="grey",
`BJ`="grey",
`GM06990`="grey",
`GM12801`="grey",
`GM12812`="grey",
`GM12813`="grey",
`GM12878`="grey",
`HUVEC`="grey",
`HeLa-S3`="grey",
`HepG2`="grey",
`IMR90`="grey",
`K562`="grey",
`MCF-7`="grey",
`NHEK`="grey",
`SK-N-SH`="grey",
`H3K27ac`='#66c2a5',
`H3K36me3`='#fc8d62',
`H3K4me1`='#8da0cb',
`H3K4me3`='#e78ac3',
`H3K9ac`='#a6d854',
`H3K27me3`='#000000',
`H3K9me3`='#e5c494'
)


make.panels <- function(tab1, tab2, colors=c(group1=1, group2=2), add.legend=FALSE, show.title=FALSE, show.xaxis=FALSE) {
print(colors)
    #if (log.y.axis == '') {
        #ylims <- range(c(tab1$enr, tab2$enr), na.rm=TRUE)*c(0.95,1.05)
        ylims <- range(pretty(c(tab1$enr, tab2$enr)))
    #} else {
    #}
print(ylims)

    # scRNAseq
    # Endothelial: there are too few endothelial cells in both scRNA and
    # scATAC-seq to have a stable signal. These should be removed.
    # Neurons: for scRNAseq only, one scRNAseq library 
    # did not have enough neurons to differentiate excitatory
    # and inhibitory because it was sorted for non-neurons. 
    plotfn(tab1[datasource == 'scrnaseq' & !(lineclass %in% c('Endothelial', 'Neurons'))],
            tab2[datasource == 'scrnaseq' & !(lineclass %in% c('Endothelial', 'Neurons'))],
            labtype='number', main='scRNA-seq', show.title=show.title,
            add.legend=add.legend, ylim=ylims, ylab='Obs / Exp',
            group1.col=colors['group1'], group2.col=colors['group2'],
            yaxis=TRUE, xaxis=show.xaxis)

    # scATAC-seq
    plotfn(tab1[datasource == 'scatacseq' & lineclass != 'Endothelial'],
            tab2[datasource == 'scatacseq' & lineclass != 'Endothelial'],
            labtype='number', main='scATAC-seq',
            add.legend=add.legend, ylim=ylims, show.title=show.title,
            group1.col=colors['group1'], group2.col=colors['group2'],
            ylab='', yaxis=FALSE, xaxis=show.xaxis)

    # RepliSeq
    plotfn(tab1[datasource == 'repliseq'],
            tab2[datasource == 'repliseq'],
            labtype='number', linetype='average',
            main='Replication timing', show.title=show.title,
            add.legend=add.legend, ylim=ylims,
            group1.col=colors['group1'], group2.col=colors['group2'],
            ylab='', yaxis=FALSE, xaxis=show.xaxis)

    # active marks
    plotfn(tab1[datasource == 'active_histone_mark'],
            tab2[datasource == 'active_histone_mark'],
            labtype='number', main='Active histone marks',
            add.legend=add.legend, ylim=ylims, show.title=show.title,
            group1.col=colors['group1'], group2.col=colors['group2'],
            ylab='', yaxis=FALSE, xaxis=show.xaxis)

    # inactive marks
    plotfn(tab1[datasource == 'inactive_histone_mark'],
            tab2[datasource == 'inactive_histone_mark'],
            labtype='number', main='Inactive histone marks',
            add.legend=add.legend, ylim=ylims, show.title=show.title,
            group1.col=colors['group1'], group2.col=colors['group2'],
            ylab='', yaxis=FALSE, xaxis=show.xaxis)
}


# if group1_col=col_map, all colors come from colmap
plotfn <- function(tab1, tab2, ylim, xlab='', ylab='', main='', typecol='lineclass', linetype=c('separate', 'average'), labtype=c('point','number'), add.legend=TRUE, col=1, yaxis=FALSE, xaxis=FALSE, group1.col=1, group2.col=2, show.title=FALSE, yaxt='s', xaxt='s', add=FALSE, use.colmap=group1.col == 'col_map', lwd=2/3, ...) {
print(group1.col)
print(group2.col)
    labtype <- match.arg(labtype)
    linetype <- match.arg(linetype)

    xlim <- range(c(tab1$quantile, tab2$quantile))
    ## changes depending on track
    tab1$type <- tab1[[typecol]]  # only allows one column for now
    tab2$type <- tab2[[typecol]]  # only allows one column for now

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

    if (use.colmap) {
        group1.lty <- lty.map[1]
        group2.lty <- lty.map[2]
    } else {
        group1.lty <- 'solid'
        group2.lty <- 'solid'
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
        #n$quantile <- as.integer(n$quantile)

        plotfn <- plot
        if (add == TRUE)
            plotfn <- lines

        # averaging over many dataclasses, so use color from the first one (arbitrarily)
        if (use.colmap) {
t <- strsplit(tab1$type[1], ' ')[[1]][1]
cat('averaging over lines, using dataclass=', t, '\n')
            group1.col <- col.map[t]
            group2.col <- col.map[t]
cat('col=', group1.col, '\n')
        }

        plotfn(tab1[,.(mean(quantile), mean(enr)), by=quantile][,.(V1,V2)],
            type='l', lwd=lwd, col=group1.col, pch=20, xlim=xlim, ylim=ylim, bty='n', xlab=xlab, ylab=ylab, main=main, xaxt=xaxt, yaxt=yaxt, lty=group1.lty, log=log.y.axis, ...)
        lines(tab2[,.(mean(quantile), mean(enr)), by=quantile][,.(V1,V2)],
            type='l', lwd=lwd, col=group2.col, pch=20, xlim=xlim, ylim=ylim, bty='n', xlab=xlab, ylab=ylab, main=main, xaxt=xaxt, yaxt=yaxt, lty=group2.lty, log=log.y.axis, ...)
        #abline(h=1, lty='dotted', col='grey', lwd=1/2)
        abline(h=1, lty='solid', col='black', lwd=1/2)
        # make a blank legend just to honor the expectation of add.legend=TRUE
        if (add.legend) {
            # separate panel because it's so large
            plot(1, pch=NA, xlab='', ylab='', xaxt='n', yaxt='n', bty='n')
        }
    } else {
        types <- sort(unique(c(tab1$type, tab2$type)))

        for (i in 1:length(types)) {
            t <- types[i]
            plot.type <- 'b'
            if (use.colmap) {
                # under col_map mode, both groups use the same color
                group1.col <- col.map[t]
                group2.col <- col.map[t]
                plot.type <- 'l'
cat("plotting trace for type t=", t, "color=", group1.col, "\n")
            }
            pf <- if (i == 1 & !add) plot else lines
            pch <- ifelse(labtype == 'point', 20, letters[i])
            pf(tab1[type == t, .(quantile, enr)], type=plot.type, lwd=lwd, col=group1.col, pch=pch,
                xlim=xlim, ylim=ylim, bty='n', xlab=xlab, ylab=ylab, main=main, xaxt=xaxt, yaxt=yaxt, lty=group1.lty, log=log.y.axis, ...)
            lines(tab2[type == t, .(quantile, enr)], type=plot.type, lwd=lwd, col=group2.col, pch=pch,
                xlim=xlim, ylim=ylim, bty='n', xlab=xlab, ylab=ylab, main=main, xaxt=xaxt, yaxt=yaxt, lty=group2.lty, log=log.y.axis, ...)
        }
        #abline(h=1, lty='dotted', col='grey', lwd=1/2)
        abline(h=1, lty='solid', col='black', lwd=1/2)

        if (add.legend) {
            # separate panel because it's so large
            par(mar=c(0,0,0,0))
            plot(1, pch=NA, xlab='', ylab='', xaxt='n', yaxt='n', bty='n')
            if (use.colmap) {
                legend('center', fill=col.map[types], legend=types, border=col.map[types], bty='n')
            } else {
                legend('center', pch=letters[1:length(types)], legend=types, bty='n')
            }
        }
    }
}


tab1 <- fread(group1.csv)
tab2 <- fread(group2.csv)

fwrite(rbind(tab1, tab2), file=out.csv)

devs=list(pdf, svglite)
outs=c(out.pdf, out.svg)
for (i in 1:2) {
    devs[[i]](width=6.5, height=(6.5/5)*3+0.1, pointsize=5, file=outs[i])

    layout(rbind(matrix(1:10, nrow=2), 10+matrix(1:10,nrow=2)), height=c(1,5,5,5))
    make.panels(tab1=tab1[muttype=='snv'], tab2=tab2[muttype=='snv'],
        colors=c(group1=group1.col, group2=group2.col),
        show.title=TRUE, show.xaxis=FALSE)
    make.panels(tab1=tab1[muttype=='indel'], tab2=tab2[muttype=='indel'],
        colors=c(group1=group1.col, group2=group2.col),
        show.title=FALSE, show.xaxis=TRUE, add.legend=TRUE)
    dev.off()
}

if ('snakemake' %in% ls()) {
    sink()
}
