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
        snakemake@output['svg']
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 6) {
    stop('usage: fig3_panel_a.R neuron_snvs.csv neuron_indels.csv oligo_snvs.csv oligo_indels.csv out.pdf out.svg')


neuron.snv <- args[1]
neuron.indel <- args[2]
oligo.snv <- args[3]
oligo.indel <- args[4]
out.pdf <- args[5]
out.svg <- args[6]

for (f in c(out.pdf, out.svg)) {
    if (file.exists(f))
        stop(paste('output file', f, 'already exists, please delete it first'))
}

suppressMessages(library(data.table))
suppressMessages(library(mutenrich))

# not using "NEUN" for neurons because "NEUN_Aug" is presumably the same
# data but more up-to-date.
n <- fread(neuron.snv)[celltype=='Excitatory-Neurons' & BINSIZE==1000 & quantile %in% 1:10 & selection != "NEUN"][order(as.integer(quantile))]
ni <- fread(neuron.indel)[celltype=='Excitatory-Neurons' & BINSIZE==1000 & quantile %in% 1:10 & selection != "NEUN"][order(as.integer(quantile))]

# Also, not using NEUN for oligos, but for a different reason: NEUN sorting
# presumably mostly removes them.
g <- fread(oligo.snv)[celltype=='Oligodendrocytes' & BINSIZE==1000 & quantile %in% 1:10 & selection != "NEUN"][order(as.integer(quantile))]
gi <- fread(oligo.indel)[celltype=='Oligodendrocytes' & BINSIZE==1000 & quantile %in% 1:10 & selection != "NEUN"][order(as.integer(quantile))]


plotfn <- function(n, g, linetype=c('separate', 'average'), labtype=c('point','number'), ...) {
    labtype <- match.arg(labtype)
    linetype <- match.arg(linetype)

    n <- n[quantile %in% 1:10]
    n$quantile <- as.integer(n$quantile)
    g <- g[quantile %in% 1:10]
    g$quantile <- as.integer(g$quantile)

    ylim <- range(n$enr, g$enr)*c(0.95,1.05)
    xlim <- range(c(n$quantile, g$quantile))
    n$type <- paste(n$donor, n$selection)
    g$type <- paste(g$donor, g$selection)

    if (linetype == 'average') {
        n$quantile <- as.integer(n$quantile)
        g$quantile <- as.integer(g$quantile)

        plot(n[,.(mean(quantile), mean(enr)), by=quantile][,.(V1,V2)],
            type='b', lwd=2, col=1, pch=20, xlim=xlim, ylim=ylim, ...)
        lines(g[,.(mean(quantile), mean(enr)), by=quantile][,.(V1,V2)],
            type='b', lwd=2, col=2, pch=20)
        abline(h=1, lty='dashed', col='grey')
    } else {
        ntypes <- unique(n$type)
        gtypes <- unique(g$type)
        for (i in 1:length(ntypes)) {
            t <- ntypes[i]
            pf <- if (i == 1) plot else lines
            pch <- ifelse(labtype == 'point', 20, letters[i])
            pf(n[type == t, .(quantile, enr)], type='b', lwd=2, col=1, pch=pch,
                xlim=xlim, ylim=ylim, ...)
        }
        for (i in 1:length(gtypes)) {
            t <- gtypes[i]
            pch <- ifelse(labtype == 'point', 20, letters[i])
            lines(g[type == t, .(quantile, enr)], type='b', lwd=2, col=2, pch=pch)
        }
        abline(h=1, lty='dashed', col='grey')

        # separate panel because it's so large
        plot(1, pch=NA, xlab='', ylab='', xaxt='n', yaxt='n', bty='n')
        legend('topleft', ncol=2, pch=c(letters[1:length(ntypes)], letters[1:length(gtypes)]),
            legend=c(ntypes,gtypes),
            col=rep(1:2, times=c(length(ntypes), length(gtypes))))
            
    }
}

layout(matrix(1:6, nrow=2, byrow=T))
plotfn(n, g, linetype='average', xlab='Expression decile', ylab='Obs/exp', main='Average enrichment SNV passA')
plotfn(n, g, labtype='number', xlab='Expression decile', ylab='Obs/exp', main='Enrichment per library SNV passA')

plotfn(ni, gi, linetype='average', xlab='Expression decile', ylab='Obs/exp', main='Average enrichment Indel passA')
plotfn(ni, gi, labtype='number', xlab='Expression decile', ylab='Obs/exp', main='Enrichment per library Indel passA')

dev.print(dev=pdf, file=out.pdf)
dev.print(dev=svg, file=out.svg)

if ('snakemake' %in% ls()) {
    sink()
}
