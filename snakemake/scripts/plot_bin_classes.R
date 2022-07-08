#!/usr/bin/env Rscript

# detect script being run by snakemake
# if so, make a mock commandArgs function
if ('snakemake' %in% ls()) {
    logfile <- snakemake@log[[1]]
    con <- file(logfile, 'w')
    sink(con, type='output')
    sink(con, type='message')

    commandArgs <- function(...) unlist(c(
        snakemake@output[1:5],  # svg/pdf output file names
        snakemake@input # variable length list of per-chromosome RDA digest files
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 5) {
    stop("usage: plot_chrom_depth.R heatmap.svg heatmap.pdf barplot.svg barplot.pdf chr_digest1.rda [ chr.digest2.rda .. chr.digestN.rda]")
}

heatmap.svg <- args[1]
heatmap.pdf <- args[2]
heatmap.jpeg <- args[3]
barplot.svg <- args[4]
barplot.pdf <- args[5]
digest.files <- args[-(1:5)]


if (file.exists(heatmap.svg))
    stop(paste('output file', heatmap.svg, 'already exists, please delete it first'))
if (file.exists(heatmap.pdf))
    stop(paste('output file', heatmap.pdf, 'already exists, please delete it first'))
if (file.exists(heatmap.jpeg))
    stop(paste('output file', heatmap.jpeg, 'already exists, please delete it first'))
if (file.exists(barplot.svg))
if (file.exists(barplot.svg))
    stop(paste('output file', barplot.svg, 'already exists, please delete it first'))
if (file.exists(barplot.pdf))
    stop(paste('output file', barplot.pdf, 'already exists, please delete it first'))

suppressMessages(library(GenomicRanges))
suppressMessages(library(GenomeInfoDb))   # For sortSeqlevels
suppressMessages(library(BSgenome))
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg19))
suppressMessages(library(extrafont))
suppressMessages(library(svglite))

if (!("Arial" %in% fonts()))
    stop("Arial font not detected; did you load extrafonts and run font_import() with the appropriate path?")

colors <- c('#df536b', '#61d04f', '#2297e6', '#28e2e5', '#cd0bbc', '#f5c710', '#9e9e9e')

gs <- lapply(digest.files, function(df) {
    cat('reading digest file', df, '\n')
    load(df)
    tiles$dp <- rowMeans(mean.mat)
    tiles.not.in.gatk$dp <- NA
    a <- c(tiles, tiles.not.in.gatk)
    a <- sort(a)
    seqlevels(a) <- paste0('chr', seqlevels(a))
    seqlevels(a) <- seqlevels(BSgenome.Hsapiens.UCSC.hg19)
    seqinfo(a) <- seqinfo(BSgenome.Hsapiens.UCSC.hg19)
    a$dpclass <- ifelse(is.na(a$dp), 0, ifelse(a$dp < 6, 1, ifelse(a$dp > quantile(a$dp, prob=0.975, na.rm=T), 2, 3)))
    a
})
names(gs) <- sapply(gs, function(g) as.character(seqnames(g)[1]))
gs <- gs[sortSeqlevels(names(gs))]


# Use a single xlim from all of the chromosomes so that they all
# plot on the same x-axis.
xlim <- range(sapply(gs, function(g) c(start(g)[1], end(g)[length(g)])))

# Get numbers of bins for each class per chromosome. For making barplot
# breakdowns of class coverage.
class.matrix <- sapply(gs, function(g) table(g$dpclass))

# Make the plots
figwidth <- 14 # inches
figheight <- 8 
# Loop over devices to save both pdf and svgs
devs <- list(svglite, pdf, jpeg)
outs <- c(heatmap.svg, heatmap.pdf, heatmap.jpeg)
for (i in 1:3) {
    if (i == 3) {
        devs[[i]](file=outs[i], width=figwidth, height=figheight, unit="in", res=600)
    } else {
        devs[[i]](file=outs[i], width=figwidth, height=figheight)
    }
    par(mar=c(5,4,1,1))
    layout(1:length(gs))
    par(mar=c(1,4,0,2), oma=c(4,0,2,0))
    for (i in 1:length(gs)) {
        g <- gs[[i]]
        image(x=start(g), y=1, z=t(t(g$dpclass)),
            bty='n',
            col=c(colors[7], colors[1], 'black', colors[2]),
            yaxt='n', ylab='', xlim=xlim,
            xaxt=ifelse(i==length(gs), 's', 'n'), xlab='')
        if (i == length(gs))
            mtext(outer=T, side=1, line=2, 'Position on chromosome')
        axis(side=2, tick=FALSE, labels=seqnames(g)[1], at=1, las=1)
    }
    dev.off()
}

plots.from.classmat <- function(class.matrix, colors) {
    rs <- rowSums(class.matrix)
    par(mar=c(5,4,3,1))
    bp <- barplot(100*rs/sum(rs), beside=T, col=colors,
        ylim=c(0,100), ylab='Percent of genome', main='Whole genome', xaxt='n')
    par(mar=c(5,2,3,1))
    text(x=bp, y=100*rs/sum(rs), pos=3,
        labels=paste0(round(100*rs/sum(rs),1), '%'))
    barplot(100*t(t(class.matrix)/colSums(class.matrix)),
        las=3, beside=T, col=colors,
        ylim=c(0,100), main='Breakdown of classes by chromosome')
}

# Make the plots
figwidth <- 14 # inches
figheight <- 8 
# Loop over devices to save both pdf and svgs
devs <- list(svglite, pdf)
outs <- c(barplot.svg, barplot.pdf)
for (i in 1:2) {
    devs[[i]](file=outs[i], width=figwidth, height=figheight)
    layout(matrix(1:4, nrow=2, byrow=T), widths=c(1,4))
    plots.from.classmat(class.matrix, colors=c(colors[7], colors[1], 'black', colors[2]))
    # row 1 is the fraction of genome absent from GATK
    plots.from.classmat(class.matrix[-1,],
        colors=c(colors[1], 'black', colors[2]))
    dev.off()
}

if ('snakemake' %in% ls()) {
    sink()
}
