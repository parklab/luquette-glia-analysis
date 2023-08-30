#!/usr/bin/env Rscript

# detect script being run by snakemake
# if so, make a mock commandArgs function
if ('snakemake' %in% ls()) {
    logfile <- snakemake@log[[1]]
    con <- file(logfile, 'w')
    sink(con, type='output')
    sink(con, type='message')

    commandArgs <- function(...) unlist(c(
        snakemake@input['sens'],
        snakemake@input['enrich'],
        snakemake@output['main_pdf'],
        snakemake@output['inset_pdf']
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}


args <- commandArgs(trailingOnly=T)

if (length(args) != 4) {
    cat("sensitivity.txt can be taken from any individual single cell; only region sizes are used from this file.\n")
    stop("usage: respfig1_genes.R sensitivity.txt enrichment.csv out.main.pdf out.inset.pdf")
}

sens.file <- args[1]
enrichment.file <- args[2]
main.pdf <- args[3]
inset.pdf <- args[4]

for (f in c(main.pdf, inset.pdf)) {
    if (file.exists(f))
        stop(sprintf("output file %s already exists, please delete it first", f))
}

# DO WORK
suppressMessages(library(scan2))
suppressMessages(library(mutenrich))


# sens.file can be from any sample, we only use the `size.mb` column,
# which records just the size of the gene tested and is independent of
# sample.
# meta==...: the sensitivity file contains enrichment tests for ALL
#     data sources that use BED regions rather than quantiles. Need to
#     restrict to just genes.
# quantile=outside: mutations not in any gene
s <- fread(sens.file)[muttype=='snv' & meta == 'datasource=gtex;dataclass=genes' & quantile != 'outside']
e <- fread(enrichment.file)[quantile != 'outside']

if (any(s$quantile != e$quantile))
    stop('sens.file and enrichment.file do not contain genes in the same order')

x <- cbind(s, e)[order(size.mb, decreasing=T)]

d <- rbindlist(lapply(seq(5000,1,-1), function(top.n) 
    x[1:top.n][, padj2 := p.adjust(pval,method='fdr')][order(padj2)][1][,.(fdr=padj2, gene=quantile, top.n=top.n)]))

# fdr < 1: when enough genes are included, all adjusted p-values become 1,
#   so the gene with the minimal FDR is a random selection from all genes.
# for marking individual genes in the inset
gd <- d[fdr < 1][order(fdr)][!duplicated(gene)]

# Main plot
pdf(file=main.pdf, width=3, height=3, pointsize=5)
plot(d$fdr, log='y', type='l', bty='n',
    ylab='-log10(min. corrected P-value)', ylim=c(0.01,1),
    xaxt='n', xlab='Number of genes considered\nordered from largest to smallest')
axis(side=1, at=pretty(1:5000), labels=rev(pretty(1:5000)))
dev.off()

# Inset
pdf(file=inset.pdf, width=2, height=2, pointsize=5)
par(mar=c(2,2,0.5,0.5))
plot(d$fdr, bty='n', xlim=c(4750,5010), log='y', type='l', xaxt='n', ylim=c(0.01,1))
axis(side=1, at=pretty(tail(1:5000,250)), labels=rev(pretty(head(1:5000,250))))
points(gd[,.(5000-top.n+1,fdr)], pch=20)
text(x=5000-gd$top.n+1, y=gd$fdr, labels=gd$gene, pos=1)
dev.off()

if ('snakemake' %in% ls()) {
    sink()
}
