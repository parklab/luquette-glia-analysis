#!/usr/bin/env Rscript

# detect script being run by snakemake
# if so, make a mock commandArgs function
if ('snakemake' %in% ls()) {
    logfile <- snakemake@log[[1]]
    con <- file(logfile, 'w')
    sink(con, type='output')
    sink(con, type='message')

    commandArgs <- function(...) unlist(c(
        snakemake@input['csv'],
        snakemake@output['pdf'],
        snakemake@output['svg'],
        snakemake@output['supppdf'],
        snakemake@output['suppsvg'],
        snakemake@output['fullsupppdf'],
        snakemake@output['fullsuppsvg']
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 7) {
    stop('usage: fig6_panel_d.R odds_ratio_matrix.csv out.pdf out.svg out.supp.pdf out.supp.svg out.fullsupp.pdf out.fullsupp.svg')
}


in.csv <- args[1]
out.pdf <- args[2]
out.svg <- args[3]
out.supp.pdf <- args[4]
out.supp.svg <- args[5]
out.fullsupp.pdf <- args[6]
out.fullsupp.svg <- args[7]

for (f in c(out.pdf, out.svg, out.supp.pdf, out.supp.svg, out.fullsupp.pdf, out.fullsupp.svg)) {
    if (file.exists(f))
        stop(paste('output file', f, 'already exists, please delete it first'))
}

suppressMessages(library(data.table))
suppressMessages(library(mutenrich))
suppressMessages(library(svglite))

odds.mat <- fread(in.csv)
types <- colnames(odds.mat)[-1]
xaxis <- odds.mat$TopNGenes

emphasize <- c('CNS-GBM','CNS-Medullo','CNS-Oligo', 'CNS-PiloAstro')
colors <- setNames(rep('grey',37), types)
colors[emphasize] <- c('orange','black','red','purple')


devs=list(pdf, svglite)
outs=c(out.pdf, out.svg)
for (i in 1:2) {
    devs[[i]](width=2, height=3.25, pointsize=5, file=outs[i])
    par(mar=c(4,8,0.1,1))
    barplot.oddsrats <- sort(setNames(as.numeric(odds.mat[TopNGenes==100,-1]), types), decreasing=FALSE)
    print(barplot.oddsrats)
    barplot(barplot.oddsrats,
        cex.names=0.8, border=F, col=colors[names(barplot.oddsrats)],
        xlim=c(0,1.4), ylab='', xlab='Odds ratio (OL odds / neuron odds)', horiz=TRUE, las=1)
    abline(v=1, lty='dashed')
    dev.off()
}


# This is the supplementary figure with a range of N (10..500) for the top N
# cancer genes.
devs=list(pdf, svglite)
outs=c(out.supp.pdf, out.supp.svg)
for (i in 1:2) {
    devs[[i]](width=6, height=3, pointsize=5, file=outs[i])
    layout(t(1:2), widths=c(3,3))  # second panel is for the massive legend

    xlim <- c(0,500)
    ylim <- c(0.8,max(as.matrix(odds.mat[,-1]), na.rm=T))
    plot(x=0,y=0, xlim=xlim, ylim=ylim, bty='n', pch=NA,
        xlab="Number of top genes",ylab="Odds ratio")
    for (i in c(2:ncol(odds.mat), which(colnames(odds.mat) %in% emphasize))) {
        col <- colors[types[i-1]]
        lwd <- ifelse(types[i-1] %in% emphasize, 2.5, 1/2)
        lines(x=xaxis, y=odds.mat[[i]], col=col, lwd=lwd)
    }

    # legend in its own panel
    plot(x=0, y=0, pch=NA, bty='n', xaxt='n', yaxt='n', xlab='', ylab='')
    legend(x="center", legend=c(emphasize, types[!(types %in% emphasize)]),
        col=c(colors[emphasize], rep('grey', length(setdiff(types, emphasize)))),
        lwd=2, bty="n", ncol=2)

    dev.off()
}

# This is the supplementary figure with a full range of N for the top N
# cancer genes.
devs=list(pdf, svglite)
outs=c(out.fullsupp.pdf, out.fullsupp.svg)
for (i in 1:2) {
    devs[[i]](width=6, height=3, pointsize=5, file=outs[i])
    layout(t(1:2), widths=c(3,3))  # second panel is for the massive legend

    xlim <- c(0,1e4)
    ylim <- c(0.5,max(as.matrix(odds.mat[,-1]), na.rm=T))
    plot(x=0,y=0, xlim=xlim, ylim=ylim, bty='n', pch=NA,
        xlab="Number of top genes",ylab="Odds ratio")
    for (i in c(2:ncol(odds.mat), which(colnames(odds.mat) %in% emphasize))) {
        col <- colors[types[i-1]]
        lwd <- ifelse(types[i-1] %in% emphasize, 2.5, 1/2)
        lines(x=xaxis, y=odds.mat[[i]], col=col, lwd=lwd)
    }

    # legend in its own panel
    plot(x=0, y=0, pch=NA, bty='n', xaxt='n', yaxt='n', xlab='', ylab='')
    legend(x="center", legend=c(emphasize, types[!(types %in% emphasize)]),
        col=c(colors[emphasize], rep('grey', length(setdiff(types, emphasize)))),
        lwd=2, bty="n", ncol=2)

    dev.off()
}

if ('snakemake' %in% ls()) {
    sink()
}
