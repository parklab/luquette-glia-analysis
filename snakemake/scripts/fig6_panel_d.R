#!/usr/bin/env Rscript

# detect script being run by snakemake
# if so, make a mock commandArgs function
if ('snakemake' %in% ls()) {
    logfile <- snakemake@log[[1]]
    con <- file(logfile, 'w')
    sink(con, type='output')
    sink(con, type='message')

    commandArgs <- function(...) unlist(c(
        snakemake@output['pdf'],
        snakemake@output['svg'],
        snakemake@output['csv'],
        snakemake@output['supppdf'],
        snakemake@output['suppsvg'],
        snakemake@output['fullsupppdf'],
        snakemake@output['fullsuppsvg'],
        snakemake@input['cancer_ors']
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 8) {
    stop('usage: fig4_panel_d.R out.pdf out.svg out.csv out.supp.pdf out.supp.svg out.fullsupp.pdf out.fullsupp.svg cancer_odds_ratios1.rda [ cancer_odd_ratios2.rda ... cancer_odds_ratiosN.rda ]')
}


out.pdf <- args[1]
out.svg <- args[2]
out.csv <- args[3]
out.supp.pdf <- args[4]
out.supp.svg <- args[5]
out.fullsupp.pdf <- args[6]
out.fullsupp.svg <- args[7]
cancer.fs <- args[-(1:7)]

for (f in c(out.pdf, out.svg, out.csv, out.supp.pdf, out.supp.svg, out.fullsupp.pdf, out.fullsupp.svg)) {
    if (file.exists(f))
        stop(paste('output file', f, 'already exists, please delete it first'))
}

suppressMessages(library(data.table))
suppressMessages(library(mutenrich))
suppressMessages(library(extrafont))
suppressMessages(library(svglite))
if (!("Arial" %in% fonts()))
    stop("Arial font not detected; did you load extrafonts and run font_import() with the appropriate path?")

data <- lapply(cancer.fs, function(f) {
    load(f)  # among other things, loads 'tops', 'types' and 'fts'
    list(tops=tops, types=types, odds=fts[1,])
})
odds.mat <- sapply(1:length(data), function(i) data[[i]]$odds)
colnames(odds.mat) <- sapply(data, function(d) d$types)
types <- colnames(odds.mat)
xaxis <- data[[1]]$tops

x <- data.table(TopNGenes=xaxis, odds.mat)[order(TopNGenes)]
fwrite(x, file=out.csv)


emphasize <- c('CNS-GBM','CNS-Medullo','CNS-Oligo', 'CNS-PiloAstro')
colors <- setNames(rep('grey',37), colnames(odds.mat))
colors[emphasize] <- c('orange','black','red','purple')


devs=list(pdf, svglite)
outs=c(out.pdf, out.svg)
for (i in 1:2) {
    devs[[i]](width=3.25, height=2, pointsize=5, file=outs[i])
    par(mar=c(8,4,3,1))
    barplot.oddsrats <- sort(odds.mat[which(xaxis==100),], decreasing=TRUE)
    print(barplot.oddsrats)
    barplot(barplot.oddsrats,
        cex.names=0.8, las=3, border=F, col=colors[names(barplot.oddsrats)],
        ylim=c(0,1.4), ylab='Odds ratio (OL odds / neuron odds)')
    abline(h=1, lty='dashed')
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
    ylim <- c(0.8,max(odds.mat, na.rm=T))
    plot(x=0,y=0, xlim=xlim, ylim=ylim, bty='n', pch=NA,
        xlab="Number of top genes",ylab="Odds ratio")
    for (i in c(1:ncol(odds.mat), which(types %in% emphasize))) {
        col <- colors[types[i]]
        lwd <- ifelse(types[i] %in% emphasize, 2.5, 1/2)
        lines(x=xaxis, y=odds.mat[,i], col=col, lwd=lwd)
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
    ylim <- c(0.5,max(odds.mat, na.rm=T))
    plot(x=0,y=0, xlim=xlim, ylim=ylim, bty='n', pch=NA,
        xlab="Number of top genes",ylab="Odds ratio")
    for (i in c(1:ncol(odds.mat), which(types %in% emphasize))) {
        col <- colors[types[i]]
        lwd <- ifelse(types[i] %in% emphasize, 2.5, 1/2)
        lines(x=xaxis, y=odds.mat[,i], col=col, lwd=lwd)
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
