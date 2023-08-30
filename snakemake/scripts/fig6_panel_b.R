#!/usr/bin/env Rscript

# detect script being run by snakemake
# if so, make a mock commandArgs function
if ('snakemake' %in% ls()) {
    logfile <- snakemake@log[[1]]
    con <- file(logfile, 'w')
    sink(con, type='output')
    sink(con, type='message')

    commandArgs <- function(...) unlist(c(
        snakemake@input['tiles'],
        snakemake@output['barplot_pdf'],
        snakemake@output['barplot_svg'],
        snakemake@output['heatmap_pdf'],
        snakemake@output['heatmap_svg'],
        snakemake@output['csv'],
        snakemake@input['cancer_qbeds'],
        snakemake@input['atac_qbeds']
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 8) {
    cat('qbed type is automatically determined by searching for datasource=scatacseq or datasource=cancer_snvdens\n')
    cat('quantiles are ignored; raw scores are used\n')
    stop('usage: fig4_panel_b.R tiles.bed barplot.pdf barplot.svg heatmap.pdf heatmap.svg out.csv qbed1 qbed2 [ qbed3 ... qbedN ]')
}


tile.file <- args[1]
barplot.pdf <- args[2]
barplot.svg <- args[3]
heatmap.pdf <- args[4]
heatmap.svg <- args[5]
out.csv <- args[6]
qbed.files <- args[-(1:6)]

for (f in c(barplot.pdf, barplot.svg, heatmap.pdf, heatmap.svg, out.csv)) {
    if (file.exists(f))
        stop(paste('output file', f, 'already exists, please delete it first'))
}

suppressMessages(library(GenomicRanges))
suppressMessages(library(pheatmap))
suppressMessages(library(data.table))
suppressMessages(library(mutenrich))
suppressMessages(library(svglite))
suppressMessages(library(grid))
suppressMessages(library(gridExtra))


gr2 <- function (bed, seqinfo = NULL, add.chr.prefix = FALSE) {
    ret <- GenomicRanges::GRanges(seqnames = bed[[1]],
        ranges = IRanges::IRanges(start = bed[[2]], bed[[3]]))
    ret$keep <- bed[,5] != 0
    ret
}

tiles <- gr2(fread(tile.file))

qbeds <- lapply(qbed.files, function(f) {
    x <- fread(f, skip=1)
    metadata <- mutenrich::read.bed.metadata(f, is.qbed=TRUE)
    if (!(metadata['datasource'] %in% c('cancer_snvdens', 'scatacseq')))
        stop('only QBEDs with datasource=cancer_snvdens or scatacseq are allowed')

    list(datasource=metadata['datasource'],
         value=ifelse(metadata['datasource'] == 'cancer_snvdens', metadata['tumor'], metadata['celltype']),
         scores=x[[5]])
})

cancer <- Filter(function(qbed) qbed$datasource == 'cancer_snvdens', qbeds)
cancer.mat <- sapply(cancer, function(x) x$scores)
colnames(cancer.mat) <- unname(sapply(cancer, function(x) x$value))

atac <- Filter(function(qbed) qbed$datasource == 'scatacseq', qbeds)
atac.mat <- sapply(atac, function(x) x$scores)
colnames(atac.mat) <- unname(sapply(atac, function(x) x$value))

str(atac.mat)
str(cancer.mat)

# Discard the lower 2.5th percentile of ATAC sites (no longer doing this, see comment below)
fit.cancer.to.atac <- function(cancer, atac, atac.eps=1e-3, cancer.eps=1e-7) {
    # IMPORTANT: multiply atac by -1 so that the 1/x xform does not invert
    # the relationship between ATAC signal and mutation density.
    data <- data.frame(atac=-1/(atac+atac.eps), cancer=(cancer+cancer.eps))
    # Removing outlier ATAC regions is no longer necessary since we are using
    # keep=TRUE tiles.
    #atac.q <- quantile(data$atac, prob=0.025, na.rm=T)
    #data <- data[tiles$keep & data$atac >= atac.q,]
    data <- data[tiles$keep,]

    lm(cancer ~ atac, data=data)
}

# Fit the model and get R squared. The linear model was useful
# for plotting/exploring, but the final value that we plot is just
# a function of the correlation.
m <- t(sapply(colnames(cancer.mat), function(cn) {
    sapply(colnames(atac.mat)[-3], function(an) {
        x <- fit.cancer.to.atac(cancer.mat[,cn], atac.mat[,an])
        summary(x)$r.squared
    })
}))

str(m)

# Get the (signed) correlation values rather than R^2
m.cor <- t(sapply(colnames(cancer.mat), function(cn) {
    sapply(colnames(atac.mat)[-3], function(an) {
        x <- fit.cancer.to.atac(cancer.mat[,cn], atac.mat[,an])
        cor(x$model$cancer, x$model$atac)
    })
}))

str(m.cor)


# Just reorder for plotting
m <- m[,c('OPC','oligo','excitatory_neuron','inhibitory_neuron','astrocyte','microglia')]
# Use prettier column names
colnames(m) <- c('OPC', 'Oligodendrocyte',
      'Excitatory-Neuron', 'Inhibitory-Neuron',
      'Astrocyte', 'Microglia')

# Just reorder for plotting
m.cor <- m.cor[,c('OPC','oligo','excitatory_neuron','inhibitory_neuron','astrocyte','microglia')]
colnames(m.cor) <- paste0(c('OPC', 'Oligodendrocyte',
      'Excitatory-Neuron', 'Inhibitory-Neuron',
      'Astrocyte', 'Microglia'), ' correlation')


# Match the colors from the UMAP plot
colors <- setNames(c('#f5c710', '#61d04f', 'black', '#28e2e5', '#cd0bbc', '#FF0000', '#9e9e9e', 'black'),
    c('Astrocyte', 'Endothelial', 'Excitatory-Neuron', 'Inhibitory-Neuron',
    'Microglia', 'Oligodendrocyte', 'OPC', 'Neuron'))


devs=list(pdf, svglite)
outs=c(barplot.pdf, barplot.svg)
for (i in 1:2) {
    devs[[i]](width=1.25, height=2, pointsize=5, file=outs[i])
    par(mar=c(8,4,3,1))
    barplot(m['CNS-GBM',], las=3, border=NA,
        ylab='R^2, GBM mutation density', col=colors[colnames(m)])
    dev.off()
}

# modified from stackoverflow:
# https://stackoverflow.com/questions/43051525/how-to-draw-pheatmap-plot-to-screen-and-also-save-to-file
save_pheatmap <- function(x, dev, filename, width=7, height=7, ...) {
   stopifnot(!missing(x))
   stopifnot(!missing(filename))
   dev(filename, width=width, height=height, ...)
   grid::grid.newpage()
   grid::grid.draw(x$gtable)
   dev.off()
}

devs=list(pdf, svglite)
outs=c(heatmap.pdf, heatmap.svg)
for (i in 1:2) {
    x <- pheatmap((m[order(m[,'OPC'],decreasing=T),]),
        main='scATACseq', cluster_row=F, cluster_cols=F, silent=TRUE, fontsize=5)
    save_pheatmap(x, dev=devs[[i]], filename=outs[i], width=2, height=4, pointsize=5)
}

# both m and m.cor are ordered by the same order(m)
fwrite(cbind(m[order(m[,'OPC'], decreasing=TRUE),], m.cor[order(m[,'OPC'], decreasing=TRUE),]), file=out.csv)

if ('snakemake' %in% ls()) {
    sink()
}
