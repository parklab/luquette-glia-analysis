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

n <- fread(neuron.snv)[eid=='E073' & BINSIZE==1000 & quantile %in% 1:10][order(as.integer(quantile))][, c('mutsfrom', 'muttype') := list('neuron', 'snv')]
ni <- fread(neuron.indel)[eid=='E073' & BINSIZE==1000 & quantile %in% 1:10][order(as.integer(quantile))][, c('mutsfrom', 'muttype') := list('neuron', 'indel')]

g <- fread(oligo.snv)[eid=='E073' & BINSIZE==1000 & quantile %in% 1:10][order(as.integer(quantile))][, c('mutsfrom', 'muttype') := list('oligo', 'snv')]
gi <- fread(oligo.indel)[eid=='E073' & BINSIZE==1000 & quantile %in% 1:10][order(as.integer(quantile))][, c('mutsfrom', 'muttype') := list('oligo', 'indel')]

outtab <- rbind(n, ni, g, gi)
fwrite(outtab, file=out.csv)

plotfn <- function(n, g, linetype=c('separate', 'average'), labtype=c('point','number'), ...) {
    labtype <- match.arg(labtype)
    linetype <- match.arg(linetype)

    n <- n[quantile %in% 1:10]
    n$quantile <- as.integer(n$quantile)
    g <- g[quantile %in% 1:10]
    g$quantile <- as.integer(g$quantile)

    ylim <- range(n$enr, g$enr)*c(0.95,1.05)
    xlim <- range(c(n$quantile, g$quantile))
    n$type <- paste(n$mark)
    g$type <- paste(g$mark)

    if (linetype == 'average') {
        n$quantile <- as.integer(n$quantile)
        g$quantile <- as.integer(g$quantile)

        plot(n[,.(mean(quantile), mean(enr)), by=quantile][,.(V1,V2)],
            type='b', lwd=2, col=1, pch=20, xlim=xlim, ylim=ylim, ...)
        lines(g[,.(mean(quantile), mean(enr)), by=quantile][,.(V1,V2)],
            type='b', lwd=2, col=2, pch=20)
        abline(h=1, lty='dashed', col='grey')
    } else {
        types <- sort(unique(c(n$type,g$type)))
        for (i in 1:length(types)) {
            t <- types[i]
            pf <- if (i == 1) plot else lines
            pch <- ifelse(labtype == 'point', 20, letters[i])
            pf(n[type == t, .(quantile, enr)], type='b', lwd=2, col=1, pch=pch,
                xlim=xlim, ylim=ylim, ...)
        }
        for (i in 1:length(types)) {
            t <- types[i]
            pch <- ifelse(labtype == 'point', 20, letters[i])
            lines(g[type == t, .(quantile, enr)], type='b', lwd=2, col=2, pch=pch)
        }
        abline(h=1, lty='dashed', col='grey')

        # separate panel because it's so large
        plot(1, pch=NA, xlab='', ylab='', xaxt='n', yaxt='n', bty='n')
        legend('topleft', pch=letters[1:length(types)], legend=types)
    }

}

devs=list(pdf, svglite)
outs=c(out.pdf, out.svg)
for (i in 1:2) {
    devs[[i]](width=5, height=3, pointsize=5, file=outs[i])
    layout(matrix(1:8, nrow=2, byrow=T))
    par(mar=c(4,4,2,1))
    plotfn(n[mark %in% c('H3K27me3', 'H3K9me3')], g[mark %in% c('H3K27me3', 'H3K9me3')], labtype='number', xlab='Pentile', ylab='Obs/exp', main='Inactive marks SNV passA\n1kb bins', family='Arial')
    plotfn(n[!(mark %in% c('H3K27me3', 'H3K9me3'))], g[!(mark %in% c('H3K27me3', 'H3K9me3'))], labtype='number', xlab='Pentile', ylab='Obs/exp', main='Active marks SNV passA\n1kb bins', family='Arial')

    plotfn(ni[mark %in% c('H3K27me3', 'H3K9me3')], gi[mark %in% c('H3K27me3', 'H3K9me3')], labtype='number', xlab='Pentile', ylab='Obs/exp', main='Inactive marks Indel passA\n1kb bins', family='Arial')
    plotfn(ni[!(mark %in% c('H3K27me3', 'H3K9me3'))], gi[!(mark %in% c('H3K27me3', 'H3K9me3'))], labtype='number', xlab='Pentile', ylab='Obs/exp', main='Active marks Indel passA\n1kb bins', family='Arial')

    dev.off()
}

if ('snakemake' %in% ls()) {
    sink()
}
