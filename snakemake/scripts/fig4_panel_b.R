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
suppressMessages(library(extrafont))
suppressMessages(library(svglite))
suppressMessages(library(grid))
suppressMessages(library(gridExtra))
if (!("Arial" %in% fonts()))
    stop("Arial font not detected; did you load extrafonts and run font_import() with the appropriate path?")


gr2 <- function (bed, seqinfo = NULL, add.chr.prefix = FALSE) {
    ret <- GenomicRanges::GRanges(seqnames = bed[[1]],
        ranges = IRanges::IRanges(start = bed[[2]], bed[[3]]))
    ret$keep <- bed[,5] != 0
    ret
}

tiles <- gr2(fread(tile.file))

qbeds <- lapply(qbed.files, function(f) {
    x <- fread(f, skip=1)
    metaline <- readLines(f, n=1)
    # 2 column matrix of k,v pairs
    metadata <- do.call(rbind, sapply(lapply(strsplit(metaline, ';')[[1]], strsplit, '='), head, 1))
    metadata <- setNames(metadata[,2], metadata[,1])
    if (!(metadata['datasource'] %in% c('cancer_snvdens', 'scatacseq')))
        stop('only QBEDs with datasource=cancer_snvdens or scatacseq are allowed')

    list(datasource=metadata['datasource'],
         value=ifelse(metadata['datasource'] == 'cancer_snvdenv', metadata['tumor'], metadata['celltype']),
         scores=x[[5]])
})

cancer <- Filter(function(qbed) qbed$datasource == 'cancer_snvdens', qbeds)
cancer.mat <- do.call(rbind, lapply(cancer, function(x) x$scores))
colnames(cancer.mat) <- sapply(cancer, function(x) x$value)

atac <- Filter(function(qbed) qbed$datasource == 'scatacseq', qbeds)
atac.mat <- do.call(rbind, lapply(atac, function(x) x$scores))
colnames(atac.mat) <- sapply(atac, function(x) x$value)

#cancer.fs <- list.files(path='/n/data1/hms/dbmi/park/jluquette/glia/analysis/try3/enrichment/cancer_snvdens/quantile/qbed',
    #pattern='cancer_snvdens___.*___normdens.1000000binsize_10quantiles.qbed',
    #full.names=T)

#atac.fs <- c(
    #'/n/data1/hms/dbmi/park/jluquette/glia/analysis/try3/enrichment/scatacseq/quantile/qbed/scatacseq___librarymerged___merged___*.1000000binsize_10quantiles.qbed',
#)


# Discard the lower 2.5th percentile of ATAC sites (=top 2.5% after 1/x xform)
fit.cancer.to.atac <- function(cancer, atac, atac.eps=1e-3, cancer.eps=1e-7) {
    data <- data.frame(atac=1/(atac+atac.eps), cancer=(cancer+cancer.eps))
    atac.q <- quantile(data$atac, prob=1-0.025, na.rm=T)
    data <- data[tiles$keep & data$atac < atac.q,]

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

# Just reorder for plotting
m <- m[,c('OPC','oligo','excitatory_neuron','inhibitory_neuron','astrocyte','microglia')]
# Use prettier column names
colnames(m) <- c('OPC', 'Oligodendrocyte',
      'Excitatory-Neuron', 'Inhibitory-Neuron',
      'Astrocyte', 'Microglia')

# Match the colors from the UMAP plot
colors <- setNames(c('#f5c710', '#61d04f', 'black', '#28e2e5', '#cd0bbc', '#FF0000', '#9e9e9e', 'black'),
    c('Astrocyte', 'Endothelial', 'Excitatory-Neuron', 'Inhibitory-Neuron',
    'Microglia', 'Oligodendrocyte', 'OPC', 'Neuron'))


devs=list(pdf, svglite)
outs=c(barplot.pdf, barplot.svg)
for (i in 1:2) {
    devs[[i]](width=2, height=2, pointsize=5, file=outs[i])
    par(mar=c(8,4,3,1))
    barplot(m['CNS-GBM',], las=3, border=NA,
        ylab='R^2, GBM mutation density', col=colors[colnames(m)])
    dev.off()
}

devs=list(pdf, svglite)
outs=c(heatmap.pdf, heatmap.svg)
for (i in 1:2) {
    devs[[i]](width=7, height=2.5, pointsize=5, file=outs[i])
    pheatmap(t(m[order(m[,'OPC'],decreasing=T),]), cluster_row=F, cluster_cols=F)
    #heatmap(t(m[order(m[,'OPC'],decreasing=T),]), Rowv=NA, Colv=NA, scale='none')
    dev.off()
}

if ('snakemake' %in% ls()) {
    sink()
}
