#!/usr/bin/env Rscript

# detect script being run by snakemake
# if so, make a mock commandArgs function
if ('snakemake' %in% ls()) {
    logfile <- snakemake@log[[1]]
    con <- file(logfile, 'w')
    sink(con, type='output')
    sink(con, type='message')

    commandArgs <- function(...) unlist(c(
        snakemake@input['neuron_muts'],
        snakemake@input['oligo_muts'],
        snakemake@input['tiles'],
        snakemake@output['pdf'],
        snakemake@output['svg'],
        snakemake@output['csv'],
        snakemake@input['cancer_bigwigs']
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 7) {
    stop('usage: fig3_panel_a.R neuron_muts.csv oligo_muts.csv genome_tiles.bed out.pdf out.svg out.csv cancer1.bigwig [ cancer2.bigwig ... cancerN.bigwig ]')
}


neuron.muts <- args[1]
oligo.muts <- args[2]
tile.file <- args[3]
out.pdf <- args[4]
out.svg <- args[5]
out.csv <- args[6]
cancer.fs <- args[-(1:6)]

for (f in c(out.pdf, out.svg, out.csv)) {
    if (file.exists(f))
        stop(paste('output file', f, 'already exists, please delete it first'))
}

suppressMessages(library(GenomicRanges))
suppressMessages(library(data.table))
suppressMessages(library(mutenrich))
suppressMessages(library(extrafont))
suppressMessages(library(svglite))
if (!("Arial" %in% fonts()))
    stop("Arial font not detected; did you load extrafonts and run font_import() with the appropriate path?")

nmut <- get(load(neuron.muts))
omut <- get(load(oligo.muts))

count.muts <- function(tiles, x)
    countOverlaps(tiles, gr(x, add.chr.prefix=TRUE))

gr2 <- function (bed, seqinfo = NULL, add.chr.prefix = FALSE) {
    ret <- GenomicRanges::GRanges(seqnames = bed[[1]],
        ranges = IRanges::IRanges(start = bed[[2]], bed[[3]]))
    ret$keep <- bed[,5] != 0
    ret$mean.dp <- bed[[6]]
    ret
}

tiles <- gr2(fread(tile.file))

nct <- count.muts(tiles, nmut)
oct <- count.muts(tiles, omut)

cancer.mat <- sapply(cancer.fs, function(f) fread(f,skip=1)[[5]])
colnames(cancer.mat) <- sapply(strsplit(colnames(cancer.mat), '___'), function(x) x[[2]])


mean.dp <- tiles$mean.dp

dp.q10 <- quantile(mean.dp[tiles$keep], prob=1/10)
dp.q90 <- quantile(mean.dp[tiles$keep], prob=9/10)


# Get correlation, R^2, p-values for cor=0 t-tests
opv <- do.call(rbind, lapply(1:ncol(cancer.mat), function(colidx) {
    col <- cancer.mat[,colidx]
    df <- data.frame(
        cancer=col[tiles$keep & mean.dp >= dp.q10 & mean.dp <= dp.q90],
        normal=oct[tiles$keep & mean.dp >= dp.q10 & mean.dp <= dp.q90])
    m <- summary(lm(cancer ~ normal, data=df))
    data.frame(Muts="Oligo",Cancer=colnames(cancer.mat)[colidx],
        Correlation=cor(df$cancer, df$normal), R.squared=m$r.squared, P.value=coef(m)['normal',4])
}))
opv <- opv[order(opv$Correlation, decreasing=TRUE),]
npv <- do.call(rbind, lapply(1:ncol(cancer.mat), function(colidx) {
    col <- cancer.mat[,colidx]
    df <- data.frame(
        cancer=col[tiles$keep & mean.dp >= dp.q10 & mean.dp <= dp.q90],
        normal=nct[tiles$keep & mean.dp >= dp.q10 & mean.dp <= dp.q90])
    m <- summary(lm(cancer ~ normal, data=df))
    data.frame(Muts="Neuron",Cancer=colnames(cancer.mat)[colidx],
        Correlation=cor(df$cancer, df$normal), R.squared=m$r.squared, P.value=coef(m)['normal',4])
}))
npv <- npv[order(npv$Correlation, decreasing=TRUE),]

fwrite(rbind(opv, npv), file=out.csv)

colors <- setNames(rep('grey',37), colnames(cancer.mat))
colors[c('CNS-GBM', 'CNS-Medullo', 'CNS-Oligo', 'CNS-PiloAstro')] <-
    c('orange','black','red','purple')

devs=list(pdf, svglite)
outs=c(out.pdf, out.svg)
for (i in 1:2) {
    devs[[i]](width=6.5, height=2, pointsize=5, file=outs[i])
    layout(t(1:2))
    par(mar=c(8,4,3,1))
    barplot(opv$Correlation, col=colors[opv$Cancer], border=F, las=3, ylim=c(-0.03,0.3),
        cex.names=0.8,
        main='Oligo passA SNVs\nno bad bins, 10% <= DP <= 90%, 1 MB bins')
    barplot(npv$Correlation, col=colors[npv$Cancer], border=F, las=3, ylim=c(-0.03,0.3),
        cex.names=0.8,
        main='Neuron passA SNVs\nno bad bins, 10% <= DP <= 90%, 1 MB bins')
    dev.off()
}

if ('snakemake' %in% ls()) {
    sink()
}
