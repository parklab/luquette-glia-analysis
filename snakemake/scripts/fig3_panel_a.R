#!/usr/bin/env Rscript

# detect script being run by snakemake
# if so, make a mock commandArgs function
if ('snakemake' %in% ls()) {
    logfile <- snakemake@log[[1]]
    con <- file(logfile, 'w')
    sink(con, type='output')
    sink(con, type='message')

    commandArgs <- function(...) unlist(c(
        snakemake@input['neuron_snv'],
        snakemake@input['neuron_indel'],
        snakemake@input['oligo_snv'],
        snakemake@input['oligo_indel'],
        snakemake@output['pdf'],
        snakemake@output['svg'],
        snakemake@output['csv']
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 7) {
    stop('usage: fig3_panel_a.R neuron_snvs.csv neuron_indels.csv oligo_snvs.csv oligo_indels.csv out.pdf out.svg out.csv')
}


neuron.snv <- args[1]
neuron.indel <- args[2]
oligo.snv <- args[3]
oligo.indel <- args[4]
out.pdf <- args[5]
out.svg <- args[6]
out.csv <- args[7]

for (f in c(out.pdf, out.svg, out.csv)) {
    if (file.exists(f))
        stop(paste('output file', f, 'already exists, please delete it first'))
}

suppressMessages(library(data.table))
suppressMessages(library(mutenrich))
suppressMessages(library(extrafont))
suppressMessages(library(svglite))
if (!("Arial" %in% fonts()))
    stop("Arial font not detected; did you load extrafonts and run font_import() with the appropriate path?")


# not using "NEUN" for neurons because "NEUN_Aug" is presumably the same
# data but more up-to-date.
# celltype==Neurons refers to only one scRNAseq library where Inh and Exc
# neurons were not separated. We don't use this library.
n <- fread(neuron.snv)[celltype != "Neurons" & BINSIZE==1000 & quantile %in% 1:100 & selection != "NEUN"][order(as.integer(quantile))][, c('mutsfrom', 'muttype') := list('neuron', 'snv')]
n$library <- paste(n$donor, n$selection)
ni <- fread(neuron.indel)[celltype != "Neurons" & BINSIZE==1000 & quantile %in% 1:100 & selection != "NEUN"][order(as.integer(quantile))][, c('mutsfrom', 'muttype') := list('neuron', 'indel')]
ni$library <- paste(ni$donor, ni$selection)

# Also, not using NEUN for oligos, but for a different reason: NEUN sorting
# presumably mostly removes them.
g <- fread(oligo.snv)[celltype != "Neurons" & BINSIZE==1000 & quantile %in% 1:100 & selection != "NEUN"][order(as.integer(quantile))][, c('mutsfrom', 'muttype') := list('oligo', 'snv')]
g$library <- paste(g$donor, g$selection)
gi <- fread(oligo.indel)[celltype != "Neurons" & BINSIZE==1000 & quantile %in% 1:100 & selection != "NEUN"][order(as.integer(quantile))][, c('mutsfrom', 'muttype') := list('oligo', 'indel')]
gi$library <- paste(gi$donor, gi$selection)

outtab <- rbind(n, ni, g, gi)
fwrite(outtab, file=out.csv)

plotfn <- function(n, g, column, linetype=c('separate', 'average'), labtype=c('point','number'), skip.legend=FALSE, ...) {
    labtype <- match.arg(labtype)
    linetype <- match.arg(linetype)

    n <- n[quantile %in% 1:10]
    n$quantile <- as.integer(n$quantile)
    g <- g[quantile %in% 1:10]
    g$quantile <- as.integer(g$quantile)

    ylim <- range(n$enr, g$enr)*c(0.95,1.05)
    xlim <- range(c(n$quantile, g$quantile))

    if (linetype == 'average') {
        n$quantile <- as.integer(n$quantile)
        g$quantile <- as.integer(g$quantile)

        plot(n[,.(mean(quantile), mean(enr)), by=quantile][,.(V1,V2)],
            type='b', lwd=2, col=1, pch=20, xlim=xlim, ylim=ylim, ...)
        lines(g[,.(mean(quantile), mean(enr)), by=quantile][,.(V1,V2)],
            type='b', lwd=2, col=2, pch=20)
        abline(h=1, lty='dashed', col='grey')
    } else {
        ntypes <- unique(n[[column]])
        gtypes <- unique(g[[column]])
        for (i in 1:length(ntypes)) {
            t <- ntypes[i]
            pf <- if (i == 1) plot else lines
            pch <- ifelse(labtype == 'point', 20, letters[i])
            ndf <- as.data.frame(n) # because using programmatic column names/values is so hard in data.table
            pf(ndf[ndf[[column]] == t, c('quantile', 'enr')], type='b', lwd=2, col=1, pch=pch,
                xlim=xlim, ylim=ylim, ...)
        }
        for (i in 1:length(gtypes)) {
            t <- gtypes[i]
            pch <- ifelse(labtype == 'point', 20, letters[i])
            gdf <- as.data.frame(g)
            lines(gdf[gdf[[column]] == t, c('quantile', 'enr')], type='b', lwd=2, col=2, pch=pch)
        }
        abline(h=1, lty='dashed', col='grey')

        if (!skip.legend) {
            # separate panel because it's so large
            plot(1, pch=NA, xlab='', ylab='', xaxt='n', yaxt='n', bty='n')
            legend('topleft', ncol=2, pch=c(letters[1:length(ntypes)], letters[1:length(gtypes)]),
                legend=c(ntypes,gtypes),
                col=rep(1:2, times=c(length(ntypes), length(gtypes))))
        }
    }
}

devs=list(pdf, svglite)
outs=c(out.pdf, out.svg)
for (i in 1:2) {
    devs[[i]](width=8, height=3, pointsize=5, file=outs[i])
    layout(matrix(1:10, nrow=2, byrow=T))
    par(mar=c(4,4,2,1))
    
    # SNVs ------------
    # For main figure plots, only using the matched cell type.
    # For average plot, don't include the combined library
    plotfn(n[donor != 'combined' & selection != 'combined' & celltype == 'Excitatory-Neurons'],
           g[donor != 'combined' & selection != 'combined' & celltype == 'Oligodendrocytes'],
        linetype='average', xlab='Expression decile', ylab='Obs/exp', main='Average enrichment SNV passAB', family='Arial')
    # For cell type specific plot, only use the combined library
    plotfn(n[donor == 'combined' & selection == 'combined'],
           g[donor == 'combined' & selection == 'combined'],
        column='celltype',
        labtype='number', xlab='Expression decile', ylab='Obs/exp', main='Enrichment per cell type SNV passAB', family='Arial', skip.legend=FALSE)
    # For per-library plot, only showing the matched cell type
    plotfn(n[celltype == 'Excitatory-Neurons'],
           g[celltype == 'Oligodendrocytes'],
        column='library',
        labtype='number', xlab='Expression decile', ylab='Obs/exp', main='Enrichment per library SNV passAB', family='Arial', skip.legend=FALSE)
    
    # Indels ------------
    plotfn(ni[donor != 'combined' & selection != 'combined' & celltype == 'Excitatory-Neurons'],
           gi[donor != 'combined' & selection != 'combined' & celltype == 'Oligodendrocytes'],
        linetype='average', xlab='Expression decile', ylab='Obs/exp', main='Average enrichment Indel passAB', family='Arial')
    plotfn(ni[donor == 'combined' & selection == 'combined'],
           gi[donor == 'combined' & selection == 'combined'],
        column='celltype',
        labtype='number', xlab='Expression decile', ylab='Obs/exp', main='Enrichment per cell type Indel passAB', family='Arial', skip.legend=FALSE)
    plotfn(ni[celltype == 'Excitatory-Neurons'],
           gi[celltype == 'Oligodendrocytes'],
        column='library',
        labtype='number', xlab='Expression decile', ylab='Obs/exp', main='Enrichment per library Indel passAB', family='Arial', skip.legend=FALSE)
    dev.off()
}

if ('snakemake' %in% ls()) {
    sink()
}
