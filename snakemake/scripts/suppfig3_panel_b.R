#!/usr/bin/env Rscript

# detect script being run by snakemake
# if so, make a mock commandArgs function
if ('snakemake' %in% ls()) {
    logfile <- snakemake@log[[1]]
    con <- file(logfile, 'w')
    sink(con, type='output')
    sink(con, type='message')

    commandArgs <- function(...) unlist(c(
        snakemake@input['neuron_repliseq'],
        snakemake@input['oligo_repliseq'],
        snakemake@output['pdf'],
        snakemake@output['svg'],
        snakemake@output['csv']
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 5) {
    stop('usage: suppfig3_panel_b.R neuron_repliseq oligo_repliseq out.pdf out.svg out.csv')
}


neuron.fs <- setNames(args[1], 'repliseq')
oligo.fs <- setNames(args[2], 'repliseq')
out.pdf <- args[3]
out.svg <- args[4]
out.csv <- args[5]

for (f in c(out.pdf, out.svg, out.csv)) {
    if (file.exists(f))
        stop(paste('output file', f, 'already exists, please delete it first'))
}

suppressMessages(library(data.table))
suppressMessages(library(mutenrich))
suppressMessages(library(svglite))


# Read in tables from various genomic covariates. These can have different
# metadata columns, so normalize them to remove metadata columns except:
#   datasource - one of: scrnaseq, scatacseq, repliseq, inactive_histone_mark, active_histone_mark
#   lineclass - for enrichment plots with multiple lines, this column defines each line
# 'ourcell' adds a column to the table recording whether these enrichment
#    analyses come from neurons or oligos.
consolidate.tables <- function(files, ourcell=c('neuron', 'oligo'), ourmuttype=c('snv','indel')) {
    ourcell <- match.arg(ourcell)
    ourmuttype <- match.arg(ourmuttype)

    wrap.fread <- function(...) {
        fread(...)[quantile %in% 1:100][order(as.integer(quantile))][, c('mutsfrom', 'muttype') := list(ourcell, ourmuttype)][BINSIZE == 1000]
    }

    # RepliSeq from ENCODE. 
    repliseq <- wrap.fread(files['repliseq'])#[celltype == 'SK-N-SH']  # shouldn't use SK-N-SH because fig4 doesn't
    # have to reverse quantile to make 1=early .. N=late
    repliseq$quantile <- (10:1)[as.integer(repliseq$quantile)]
    repliseq$lineclass <- paste(repliseq$celltype, repliseq$geoid)
    repliseq$celltype <- NULL
    repliseq$geoid <- NULL

    repliseq
}


make.panels <- function(tab, sig, colors=c(neuron=1, oligo=2), add.legend=FALSE, show.title=FALSE, show.xaxis=FALSE) {
    ylims <- range(tab$enr)*c(0.95,1.05)
print(ylims)

    plotfn(tab[datasource == 'repliseq'],
            labtype='point', linetype='separate',
            main=sprintf('ENCODE RepliSeq\n(%d tissues)', length(unique(tab$lineclass))),
            show.title=TRUE,
            add.legend=TRUE, ylim=ylims,
            ylab='Obs/Exp', xlab='Replication timing decile',
            xaxis=TRUE, yaxis=TRUE)
}



plotfn <- function(n, ylim, xlab='', ylab='', main='', typecol='lineclass', linetype=c('separate', 'average'), labtype=c('point','number'), add.legend=TRUE, col=1, yaxis=FALSE, xaxis=FALSE, neuron.col=1, oligo.col=2, show.title=FALSE, yaxt='s', xaxt='s', add=FALSE, ...) {
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

        plotfn(n[mutsfrom == 'neuron',.(mean(quantile), mean(enr)), by=quantile][,.(V1,V2)],
            type='b', lwd=1/2, col=neuron.col, pch=20, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, main=main, xaxt=xaxt, yaxt=yaxt, ...)
        lines(n[mutsfrom == 'oligo',.(mean(quantile), mean(enr)), by=quantile][,.(V1,V2)],
            type='b', lwd=1/2, col=oligo.col, pch=20, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, main=main, xaxt=xaxt, yaxt=yaxt, ...)
        abline(h=1, lty='dashed', col='grey')
        # make a blank legend just to honor the expectation of add.legend=TRUE
        if (add.legend) {
            # separate panel because it's so large
            plot(1, pch=NA, xlab='', ylab='', xaxt='n', yaxt='n', bty='n')
        }
    } else {
        types <- sort(unique(n$type))

        for (i in 1:length(types)) {
            t <- types[i]
            pf <- if (i == 1 & !add) plot else lines
            pch <- ifelse(labtype == 'point', 20, letters[i])
            pf(n[mutsfrom == 'neuron' & type == t, .(quantile, enr)], type='b', lwd=1/2, col=neuron.col, pch=pch,
                xlim=xlim, ylim=ylim, bty='n', xlab=xlab, ylab=ylab, main=main, xaxt=xaxt, yaxt=yaxt, ...)
            lines(n[mutsfrom == 'oligo' & type == t, .(quantile, enr)], type='b', lwd=1/2, col=oligo.col, pch=pch,
                xlim=xlim, ylim=ylim, bty='n', xlab=xlab, ylab=ylab, main=main, xaxt=xaxt, yaxt=yaxt, ...)
        }
        abline(h=1, lty='dashed', col='grey')

        if (add.legend) {
            # separate panel because it's so large
            par(mar=c(0,0,0,0))
            plot(1, pch=NA, xlab='', ylab='', xaxt='n', yaxt='n', bty='n')
            legend('topleft', pch=letters[1:length(types)], legend=types, bty='n')
        }
    }
}


n <- consolidate.tables(files=neuron.fs, ourcell='neuron')
o <- consolidate.tables(files=oligo.fs, ourcell='oligo')
all.data <- rbind(n, o)

fwrite(all.data, file=out.csv)

devs=list(pdf, svglite)
outs=c(out.pdf, out.svg)
for (i in 1:2) {
    devs[[i]](width=1.5, height=2.5, pointsize=5, file=outs[i])
    layout(1:3, heights=c(1/3,4,2.75))
    make.panels(tab=all.data, show.title=TRUE)
    dev.off()
}

if ('snakemake' %in% ls()) {
    sink()
}
