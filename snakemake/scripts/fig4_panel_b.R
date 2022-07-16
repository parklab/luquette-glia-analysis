#!/usr/bin/env Rscript

library(GenomicRanges)
library(mutenrich)
library(data.table)
library(svglite)
library(pheatmap)


gr2 <- function (bed, seqinfo = NULL, add.chr.prefix = FALSE) {
    ret <- GenomicRanges::GRanges(seqnames = bed[[1]],
        ranges = IRanges::IRanges(start = bed[[2]], bed[[3]]))
    ret$keep <- bed[,5] != 0
    ret
}

tiles <- gr2(fread('/n/data1/hms/dbmi/park/jluquette/glia/analysis/try3/alignability/genome_tiles/genome_tiles_1000000binsize.bed'))

cancer.fs <- list.files(path='/n/data1/hms/dbmi/park/jluquette/glia/analysis/try3/enrichment/cancer_snvdens/quantile/qbed',
    pattern='cancer_snvdens___.*___normdens.1000000binsize_10quantiles.qbed',
    full.names=T)

cancer.mat <- sapply(cancer.fs, function(f) fread(f,skip=1)[[5]])
colnames(cancer.mat) <- sapply(strsplit(colnames(cancer.mat), '___'), function(x) x[[2]])



atac.fs <- c(
    '/n/data1/hms/dbmi/park/jluquette/glia/analysis/try3/enrichment/scatacseq/quantile/qbed/scatacseq___librarymerged___merged___OPC.1000000binsize_10quantiles.qbed',
    '/n/data1/hms/dbmi/park/jluquette/glia/analysis/try3/enrichment/scatacseq/quantile/qbed/scatacseq___librarymerged___merged___astrocyte.1000000binsize_10quantiles.qbed',
    '/n/data1/hms/dbmi/park/jluquette/glia/analysis/try3/enrichment/scatacseq/quantile/qbed/scatacseq___librarymerged___merged___endothelial.1000000binsize_10quantiles.qbed',
    '/n/data1/hms/dbmi/park/jluquette/glia/analysis/try3/enrichment/scatacseq/quantile/qbed/scatacseq___librarymerged___merged___excitatory_neuron.1000000binsize_10quantiles.qbed',
    '/n/data1/hms/dbmi/park/jluquette/glia/analysis/try3/enrichment/scatacseq/quantile/qbed/scatacseq___librarymerged___merged___inhibitory_neuron.1000000binsize_10quantiles.qbed',
    '/n/data1/hms/dbmi/park/jluquette/glia/analysis/try3/enrichment/scatacseq/quantile/qbed/scatacseq___librarymerged___merged___microglia.1000000binsize_10quantiles.qbed',
    '/n/data1/hms/dbmi/park/jluquette/glia/analysis/try3/enrichment/scatacseq/quantile/qbed/scatacseq___librarymerged___merged___oligo.1000000binsize_10quantiles.qbed'
)
atac.fc.fs <- c(
    '/n/data1/hms/dbmi/park/jluquette/glia/analysis/negbin_model/foldchange_qbedenrich/scatacseq_OPC.1000000binsize_10quantiles.qbed',
    '/n/data1/hms/dbmi/park/jluquette/glia/analysis/negbin_model/foldchange_qbedenrich/scatacseq_astrocyte.1000000binsize_10quantiles.qbed',
    '/n/data1/hms/dbmi/park/jluquette/glia/analysis/negbin_model/foldchange_qbedenrich/scatacseq_excitatory_neuron.1000000binsize_10quantiles.qbed',
    '/n/data1/hms/dbmi/park/jluquette/glia/analysis/negbin_model/foldchange_qbedenrich/scatacseq_inhibitory_neuron.1000000binsize_10quantiles.qbed',
    '/n/data1/hms/dbmi/park/jluquette/glia/analysis/negbin_model/foldchange_qbedenrich/scatacseq_microglia.1000000binsize_10quantiles.qbed',
    '/n/data1/hms/dbmi/park/jluquette/glia/analysis/negbin_model/foldchange_qbedenrich/scatacseq_oligo.1000000binsize_10quantiles.qbed')


# For absolute atac signal
atac.mat <- sapply(atac.fs, function(f) fread(f,skip=1)[[5]])
colnames(atac.mat) <- sapply(strsplit(sapply(strsplit(colnames(atac.mat), '___'), `[`, 4), '.', fixed=TRUE), head, 1)

# For foldchange atac signal
atac.fc.mat <- sapply(atac.fc.fs, function(f) fread(f,skip=1)[[5]])
colnames(atac.fc.mat) <- sub(pattern='scatacseq_', replacement='', x=sapply(strsplit(basename(atac.fc.fs), '.', fixed=TRUE), head, 1))



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

# I sized this manually via X11
pheatmap(m[order(m[,'OPC'],decreasing=T),], cluster_row=F, cluster_cols=F)
#dev.print(dev=pdf, file='cancer_mut_predicted_by_atac.1mb_nobadbins_atac_gt_2.5pctile.pdf')
#dev.print(dev=svglite, file='cancer_mut_predicted_by_atac.1mb_nobadbins_atac_gt_2.5pctile.svglite')

colors <- setNames(c('#f5c710', '#61d04f', 'black', '#28e2e5', '#cd0bbc', '#FF0000', '#9e9e9e', 'black'),
    c('Astrocyte', 'Endothelial', 'Excitatory-Neuron', 'Inhibitory-Neuron',
      'Microglia', 'Oligodendrocyte', 'OPC', 'Neuron'))


# Also manually sized with X11
par(mar=c(8,4,3,1))
barplot(m['CNS-GBM',], las=3, border=NA,
    ylab='R^2, GBM mutation density', col=colors[colnames(m)])
#dev.print(dev=pdf, file='cancer_mut_predicted_by_atac.1mb_nobadbins_atac_gt_2.5pctile.JUST_GBM.pdf')
#dev.print(dev=svglite, file='cancer_mut_predicted_by_atac.1mb_nobadbins_atac_gt_2.5pctile.JUST_GBM.svglite')


# Wasn't too interesting.
if (FALSE) {
# See what foldchange looks like
mfc <- t(sapply(colnames(cancer.mat), function(cn) {
    sapply(colnames(atac.fc.mat), function(an) {
        x <- fit.cancer.to.atac(cancer.mat[,cn], atac.fc.mat[,an])
        summary(x)$r.squared
    })
}))

# Just reorder for plotting
mfc <- mfc[,c('OPC','oligo','excitatory_neuron','inhibitory_neuron','astrocyte','microglia')]
# Use prettier column names
colnames(mfc) <- c('OPC', 'Oligodendrocyte',
      'Excitatory-Neuron', 'Inhibitory-Neuron',
      'Astrocyte', 'Microglia')

# I sized this manually via X11
pheatmap(mfc[order(mfc[,'OPC'],decreasing=T),], cluster_row=F, cluster_cols=F)
dev.print(dev=pdf, file='cancer_mut_predicted_by_atac_FOLD_CHANGE.1mb_nobadbins_atac_gt_2.5pctile.pdf')
dev.print(dev=svglite, file='cancer_mut_predicted_by_atac_FOLD_CHANGE.1mb_nobadbins_atac_gt_2.5pctile.svglite')
}
