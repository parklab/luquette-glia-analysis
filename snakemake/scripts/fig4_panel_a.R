#!/usr/bin/env Rscript

library(GenomicRanges)
library(mutenrich)
library(data.table)
library(svglite)

nmut <- get(load('/n/data1/hms/dbmi/park/jluquette/glia/analysis/try3/input/neuron___mut___A.rda'))
omut <- get(load('/n/data1/hms/dbmi/park/jluquette/glia/analysis/try3/input/oligo___mut___A.rda'))


count.muts <- function(tiles, x)
    countOverlaps(tiles, gr(x, add.chr.prefix=TRUE))

gr2 <- function (bed, seqinfo = NULL, add.chr.prefix = FALSE) {
    ret <- GenomicRanges::GRanges(seqnames = bed[[1]],
        ranges = IRanges::IRanges(start = bed[[2]], bed[[3]]))
    ret$keep <- bed[,5] != 0
    ret
}

tiles <- gr2(fread('/n/data1/hms/dbmi/park/jluquette/glia/analysis/try3/alignability/genome_tiles/genome_tiles_1000000binsize.bed'))

nct <- count.muts(tiles, nmut)
oct <- count.muts(tiles, omut)

cancer.fs <- list.files(path='/n/data1/hms/dbmi/park/jluquette/glia/analysis/try3/enrichment/cancer_snvdens/quantile/qbed',
    pattern='cancer_snvdens___.*___normdens.1000000binsize_10quantiles.qbed',
    full.names=T)

cancer.mat <- sapply(cancer.fs, function(f) fread(f,skip=1)[[5]])
colnames(cancer.mat) <- sapply(strsplit(colnames(cancer.mat), '___'), function(x) x[[2]])


# Get mean depth profile across single cells
depth.fs <- list.files(path='/n/data1/hms/dbmi/park/jluquette/glia/analysis/try3/enrichment/depth/quantile/qbed/',
    pattern='depth___.*___.*.1000000binsize_50quantiles.qbed', full.names=T)

# instead of sensitivity, try technical covariates as proxies
depth <- do.call(data.frame, lapply(depth.fs, function(f) fread(f, skip=1)[[5]]))
colnames(depth) <- sub('.1000000binsize_50quantiles.qbed', '',
    sapply(strsplit(depth.fs, split='___', fixed=TRUE), function(x) x[3]))

mean.dp <- rowMeans(depth)

dp.q10 <- quantile(mean.dp[tiles$keep], prob=1/10)
dp.q90 <- quantile(mean.dp[tiles$keep], prob=9/10)

ocor <- sort(apply(cancer.mat, 2, function(col)
    cor(col[tiles$keep & mean.dp >= dp.q10 & mean.dp <= dp.q90],
        oct[tiles$keep & mean.dp >= dp.q10 & mean.dp <= dp.q90])))
ncor <- sort(apply(cancer.mat, 2, function(col)
    cor(col[tiles$keep & mean.dp >= dp.q10 & mean.dp <= dp.q90],
        nct[tiles$keep & mean.dp >= dp.q10 & mean.dp <= dp.q90])))


colors <- setNames(rep('grey',37), colnames(cancer.mat))
colors[c('CNS-GBM', 'CNS-Medullo', 'CNS-Oligo', 'CNS-PiloAstro')] <-
    c('orange','black','red','purple')


layout(1:2)
par(mar=c(8,4,3,1))
barplot(ocor, col=colors[names(ocor)], border=F, las=3, ylim=c(-0.03,0.3),
    cex.names=0.8,
    main='Oligo passA SNVs\nno bad bins, 10% <= DP <= 90%, 1 MB bins')
barplot(ncor, col=colors[names(ncor)], border=F, las=3, ylim=c(-0.03,0.3),
    cex.names=0.8,
    main='Neuron passA SNVs\nno bad bins, 10% <= DP <= 90%, 1 MB bins')

# NB. I made the size of these plots manually
#dev.print(dev=pdf, file='cancer_correlation_normdens_1000000binsize_tileskeep_dp_gt_10pct_and_lt_90pct.pdf')
#dev.print(dev=svglite, file='cancer_correlation_normdens_1000000binsize_tileskeep_dp_gt_10pct_and_lt_90pct.svg')


# Get p-values for cor=0 t-tests
opv <- apply(cancer.mat, 2, function(col) {
    m <- summary(lm(cancer ~ normal, data=data.frame(
        cancer=col[tiles$keep & mean.dp >= dp.q10 & mean.dp <= dp.q90],
        normal=oct[tiles$keep & mean.dp >= dp.q10 & mean.dp <= dp.q90])))
    c(m$r.squared, coef(m)['normal',4])
})
npv <- apply(cancer.mat, 2, function(col) {
    m <- summary(lm(cancer ~ normal, data=data.frame(
        cancer=col[tiles$keep & mean.dp >= dp.q10 & mean.dp <= dp.q90],
        normal=nct[tiles$keep & mean.dp >= dp.q10 & mean.dp <= dp.q90])))
    c(m$r.squared, coef(m)['normal',4])
})
