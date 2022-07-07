#!/usr/bin/env Rscript

# detect script being run by snakemake
# if so, make a mock commandArgs function
if ('snakemake' %in% ls()) {
    logfile <- snakemake@log[[1]]
    con <- file(logfile, 'w')
    sink(con, type='output')
    sink(con, type='message')

    commandArgs <- function(...) unlist(c(
        snakemake@input[1], snakemake@output[1:3]
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 4) {
    stop("usage: plot_atac_umap.r seurat_umap.csv out.svg out.pdf out.jpeg")
}

in.csv <- args[1]
out.svg <- args[2]
out.pdf <- args[3]
out.jpeg <- args[4]

if (file.exists(out.svg))
    stop(paste('output file', out.svg, 'already exists, please delete it first'))
if (file.exists(out.pdf))
    stop(paste('output file', out.pdf, 'already exists, please delete it first'))
if (file.exists(out.jpeg))
    stop(paste('output file', out.jpeg, 'already exists, please delete it first'))

suppressMessages(library(extrafont))
suppressMessages(library(svglite))
suppressMessages(library(data.table))

if (!("Arial" %in% fonts()))
    stop("Arial font not detected; did you load extrafonts and run font_import() with the appropriate path?")

cat("Reading", in.csv, "\n")
d <- fread(in.csv)

figwidth <- 6.25
figheight <- 4.5

# uses a yellow color that's hard to see
#colors <- 1:7
#colors <- c('#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#a65628','#f781bf')
#colors <- c('#df536b', '#61d04f', '#2297e6', '#28e2e5', '#cd0bbc', '#f5c710', '#9e9e9e')
# now to match neuron=black, oligo=red
colors <- c('#f5c710', '#61d04f', 'black', '#28e2e5', '#cd0bbc', '#FF0000', '#9e9e9e')
# case insensitive sorting
names(colors) <- names(table(d$CellType))[order(tolower(names(table(d$CellType))))]

# Loop over devices to save both pdf and svgs
devs <- list(svglite, pdf, jpeg)
outs <- c(out.svg, out.pdf, out.jpeg)
for (i in 1:3) {
    if (i == 3) {
        devs[[i]](file=outs[i], width=figwidth, height=figheight, unit="in", res=600)
    } else {
        devs[[i]](file=outs[i], width=figwidth, height=figheight)
    }

    layout(t(1:2), width=c(2.75,1))
    par(mar=c(4,4,1,0.1))
    plot(d$UMAP1, d$UMAP2, pch='.', col=colors[d$CellType],
        xlab='UMAP 1', ylab='UMAP 2')
    # hack to get the legend plotted in a separate panel
    par(mar=c(0.1,0.1,0.1,0.1))
    plot(1, xaxt='n', xlab='', yaxt='n', ylab='', main='', bty='n', pch=NA)
    legend('center', fill=colors, legend=names(colors), bty='n')

    dev.off()
}


if ('snakemake' %in% ls()) {
    sink()
}
