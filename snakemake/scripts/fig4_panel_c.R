suppressMessages(library(mutenrich))
suppressMessages(library(extrafont))
suppressMessages(library(svglite))


if (!("Arial" %in% fonts()))
    stop("Arial font not detected; did you load extrafonts and run font_import() with the appropriate path?")


# 1MB plots
#load('scatacseq_foldchange_oligo_AB_1000000.SUMMARY.rda'); aoes.1m <- lapply(es, function(e) e[as.character(1:10),])
#load('scatacseq_foldchange_neuron_AB_1000000.SUMMARY.rda'); anes.1m <- lapply(es, function(e) e[as.character(1:10),])
load('scatacseq_foldchange_oligo_AB_1000.SUMMARY.rda'); aoes.1m <- lapply(es, function(e) e[as.character(1:10),])
load('scatacseq_foldchange_neuron_AB_1000.SUMMARY.rda'); anes.1m <- lapply(es, function(e) e[as.character(1:10),])
atacmeta <- emeta
atacmap <- setNames(c('OPCs', 'Astrocytes', 'Excitatory-Neurons', 'Inhibitory-Neurons',
    'Microglia', 'Oligodendrocytes'),
    c('OPC', 'astrocyte', 'excitatory_neuron', 'inhibitory_neuron', 'microglia', 'oligo'))
print(str(atacmeta$celltype))
print(str(levels(atacmeta$celltype)))
print(atacmeta$celltype)
atacmeta$celltype <- atacmap[levels(atacmeta$celltype)]
print(atacmeta$celltype)

# 1 MB bins
#load('scrnaseq_foldchange_oligo_AB_1000000.SUMMARY.rda'); roes.1m <- lapply(es, function(e) e[as.character(1:10),])
#load('scrnaseq_foldchange_neuron_AB_1000000.SUMMARY.rda'); rnes.1m <- lapply(es, function(e) e[as.character(1:10),])
load('scrnaseq_foldchange_oligo_AB_1000.SUMMARY.rda'); roes.1m <- lapply(es, function(e) e[as.character(1:10),])
load('scrnaseq_foldchange_neuron_AB_1000.SUMMARY.rda'); rnes.1m <- lapply(es, function(e) e[as.character(1:10),])
rnameta <- emeta
rnameta$celltype <- levels(rnameta$celltype)


# Matches the UMAP figure colors
colors <- setNames(c('#f5c710', '#61d04f', 'black', '#28e2e5', '#cd0bbc', '#FF0000', '#9e9e9e', 'black'),
    c('Astrocytes', 'Endothelial', 'Excitatory-Neurons', 'Inhibitory-Neurons',
      'Microglia', 'Oligodendrocytes', 'OPCs', 'Neurons'))

textpane <- function(txt, ...) {
    oldmar <- par(mar=c(0,0,0,0))$mar
    plot(1, pch=NA, bty='n', xaxt='n', yaxt='n')
    legend('center', legend=txt, bty='n', cex=1.4, ...)
    par(mar=oldmar)
}


x11(width=9.5, height=8)
layout(matrix(1:6, nrow=2,byrow=T), width=c(2,2,1))
n=5
yl=c(0.5,1.5)

# For snRNAseq, exc. and inh. neurons were merged. Can break this
# out later.
for (i in 1:n)
    plot.enrich(es=rnes.1m[[i]], main='Neuron SNVs passAB 1MB bins', type='l', ltype='l',
        xlab='snRNAseq foldchange (decile)', las=2, 
        lcol=colors[rnameta$celltype[i]], add=i>1, ylim=yl, bootstrap.ci=F)
abline(h=1, lty='dashed', col='grey')
for (i in 1:n)
    plot.enrich(es=roes.1m[[i]], main='Oligodendrocyte SNVs passAB 1MB bins', type='l', ltype='l',
        xlab='snRNAseq foldchange (decile)', las=2,
        lcol=colors[rnameta$celltype[i]], add=i>1, ylim=yl, bootstrap.ci=F)
abline(h=1, lty='dashed', col='grey')

textpane(x='center', txt=rnameta$celltype, lwd=2, col=colors[rnameta$celltype])

#dev.print(dev=pdf, file='scrnaseq_foldchange_plots_1m_passAB.pdf')
#dev.print(dev=svglite, file='scrnaseq_foldchange_plots_1m_passAB.svg')

#x11(width=9.5, height=4)
#layout(t(1:3), width=c(2,2,1))
n=6
# for 1 MB plots
#yl=c(0.5,1.5)
yl=c(0.85,1.2)  # for 1 KB plots

# For snRNAseq, exc. and inh. neurons were merged. Can break this
# out later.
for (i in 1:n)
    plot.enrich(es=anes.1m[[i]], main='Neuron SNVs passAB 1MB bins', type='l', ltype='l',
        xlab='scATACseq foldchange (decile)', las=2, 
        lcol=colors[atacmeta$celltype[i]], add=i>1, ylim=yl, bootstrap.ci=F)
abline(h=1, lty='dashed', col='grey')
for (i in 1:n)
    plot.enrich(es=aoes.1m[[i]], main='Oligodendrocyte SNVs passAB 1MB bins', type='l', ltype='l',
        xlab='scATACseq foldchange (decile)', las=2,
        lcol=colors[atacmeta$celltype[i]], add=i>1, ylim=yl, bootstrap.ci=F)
abline(h=1, lty='dashed', col='grey')

textpane(x='center', txt=atacmeta$celltype, lwd=2, col=colors[atacmeta$celltype])

#dev.print(dev=pdf, file='scatacseq_foldchange_plots_1m_passAB.pdf')
#dev.print(dev=svglite, file='scatacseq_foldchange_plots_1m_passAB.svg')

