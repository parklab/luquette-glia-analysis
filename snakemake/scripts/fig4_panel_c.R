#!/usr/bin/env Rscript

# detect script being run by snakemake
# if so, make a mock commandArgs function
if ('snakemake' %in% ls()) {
    logfile <- snakemake@log[[1]]
    con <- file(logfile, 'w')
    sink(con, type='output')
    sink(con, type='message')

    commandArgs <- function(...) unlist(c(
        snakemake@input['neuron_rna'],
        snakemake@input['neuron_atac'],
        snakemake@input['oligo_rna'],
        snakemake@input['oligo_atac'],
        snakemake@output['pdf'],
        snakemake@output['svg'],
        snakemake@output['csv']
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 7) {
    stop('usage: fig4_panel_c.R neuron_rna.rda neuron_atac.rda oligo_rna.rda oligo_atac.rda out.pdf out.svg out.csv')
}

neuron.rna <- args[1]
neuron.atac <- args[2]
oligo.rna <- args[3]
oligo.atac <- args[4]
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

# ATACseq
load(oligo.atac)
oatactab <- do.call(rbind, lapply(1:nrow(emeta), function(i) { x <- data.frame(Assay='scATACseq', MutsFrom='Oligo', AssayCelltype=emeta$celltype[i], Quantile=rownames(es[[i]]), es[[i]]); x[as.character(1:10),] }))
aoes <- lapply(es, function(e) e[as.character(1:10),])
load(neuron.atac)
natactab <- do.call(rbind, lapply(1:nrow(emeta), function(i) { x <- data.frame(Assay='scATACseq', MutsFrom='Neuron', AssayCelltype=emeta$celltype[i], Quantile=rownames(es[[i]]), es[[i]]); x[as.character(1:10),] }))
anes <- lapply(es, function(e) e[as.character(1:10),])
atacmeta <- emeta
atacmap <- setNames(c('OPCs', 'Astrocytes', 'Excitatory-Neurons', 'Inhibitory-Neurons',
    'Microglia', 'Oligodendrocytes', 'Endothelial'),
    c('OPC', 'astrocyte', 'excitatory_neuron', 'inhibitory_neuron', 'microglia', 'oligo', 'endothelial'))
print(atacmeta$celltype)
atacmeta$celltype <- atacmap[atacmeta$celltype]
print(atacmeta$celltype)

# RNAseq
load(oligo.rna)
ornatab <- do.call(rbind, lapply(1:nrow(emeta), function(i) { x <- data.frame(Assay='scRNAseq', MutsFrom='Oligo', AssayCelltype=emeta$celltype[i], Quantile=rownames(es[[i]]), es[[i]]); x[as.character(1:10),] }))
roes <- lapply(es, function(e) e[as.character(1:10),])
load(neuron.rna)
nrnatab <- do.call(rbind, lapply(1:nrow(emeta), function(i) { x <- data.frame(Assay='scRNAseq', MutsFrom='Neuron', AssayCelltype=emeta$celltype[i], Quantile=rownames(es[[i]]), es[[i]]); x[as.character(1:10),] }))
rnes <- lapply(es, function(e) e[as.character(1:10),])
rnameta <- emeta


fwrite(rbind(oatactab,ornatab, natactab,nrnatab), file=out.csv)

# Match the UMAP figure colors
colors <- setNames(c('#f5c710', '#61d04f', 'black', '#28e2e5', '#cd0bbc', '#FF0000', '#9e9e9e', 'black'),
    c('Astrocytes', 'Endothelial', 'Excitatory-Neurons', 'Inhibitory-Neurons',
      'Microglia', 'Oligodendrocytes', 'OPCs', 'Neurons'))

textpane <- function(txt, ...) {
    oldmar <- par(mar=c(0,0,0,0))$mar
    plot(1, pch=NA, bty='n', xaxt='n', yaxt='n')
    legend('center', legend=txt, bty='n', cex=1.4, ...)
    par(mar=oldmar)
}


devs=list(pdf, svglite)
outs=c(out.pdf, out.svg)
for (i in 1:2) {
    devs[[i]](width=4.5, height=2.5, pointsize=5, file=outs[i])

    # format:
    # row1: scRNAseq: neuron plot, oligo plot, legend 
    # row2: scATACseq: neuron plot, oligo plot, legend 
    layout(matrix(1:6, nrow=2,byrow=T), width=c(2,2,1))
    par(mar=c(4,4,2,1))

    n=6
    yl=c(0.5,1.5)

    # For snRNAseq, exc. and inh. neurons were merged. Can break this
    # out in supplement if we want.
    for (i in 1:n)
        plot.enrich(es=rnes[[i]], main='Neuron SNVs', type='l', ltype='b',
            xlab='snRNAseq foldchange (decile)',
            lcol=colors[rnameta$celltype[i]], add=i>1, ylim=yl, bootstrap.ci=F)
    abline(h=1, lty='dashed', col='grey')
    for (i in 1:n)
        plot.enrich(es=roes[[i]], main='Oligodendrocyte SNVs', type='l', ltype='b',
            xlab='snRNAseq foldchange (decile)',
            lcol=colors[rnameta$celltype[i]], add=i>1, ylim=yl, bootstrap.ci=F)
    abline(h=1, lty='dashed', col='grey')

    textpane(x='center', txt=rnameta$celltype, lwd=2, col=colors[rnameta$celltype])

    n=7
    yl=c(0.5,1.5) # for 1 MB plots
    #yl=c(0.85,1.2)  # for 1 KB plots

    for (i in 1:n)
        plot.enrich(es=anes[[i]], main='Neuron SNVs', type='l', ltype='b',
            xlab='scATACseq foldchange (decile)',
            lcol=colors[atacmeta$celltype[i]], add=i>1, ylim=yl, bootstrap.ci=F)
    abline(h=1, lty='dashed', col='grey')
    for (i in 1:n)
        plot.enrich(es=aoes[[i]], main='Oligodendrocyte SNVs', type='l', ltype='b',
            xlab='scATACseq foldchange (decile)',
            lcol=colors[atacmeta$celltype[i]], add=i>1, ylim=yl, bootstrap.ci=F)
    abline(h=1, lty='dashed', col='grey')

    textpane(x='center', txt=atacmeta$celltype, lwd=2, col=colors[atacmeta$celltype])

    dev.off()
}

if ('snakemake' %in% ls()) {
    sink()
}
