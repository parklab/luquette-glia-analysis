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
        snakemake@output['inset_pdf'],
        snakemake@output['all_csv'],
        snakemake@output['minpoints_csv']
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}


args <- commandArgs(trailingOnly=T)

if (length(args) != 6) {
    cat("sensitivity.txt can be taken from any individual single cell; only region sizes are used from this file.\n")
    stop("usage: respfig1_genes.R sensitivity.txt enrichment.csv out.main.pdf out.inset.pdf")
}

sens.file <- args[1]
enrichment.file <- args[2]
main.pdf <- args[3]
inset.pdf <- args[4]
all.csv <- args[5]
minpoints.csv <- args[6]

for (f in c(main.pdf, inset.pdf, all.csv, minpoints.csv)) {
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
s <- fread(sens.file)[metaline == 'datasource=gtex;dataclass=genes' & feature != 'outside']
e <- fread(enrichment.file)[quantile != 'outside']

x <- e[s,,on=.(quantile=feature)][order(width, decreasing=TRUE)]

n.genes <- nrow(x)
#n.genes <- 1300      # small test set
d <- rbindlist(lapply(n.genes:1, function(top.n) 
    x[1:top.n][, padj2 := p.adjust(pval,method='fdr')][order(padj2)][1][,.(fdr=padj2, gene=quantile, top.n=top.n)]))

fwrite(d, file=all.csv)

# gd gets each gene that attains the minimum P-value and the earliest such cutoff
# N at which the gene attains the minimum.  This is useful for plotting the
# individual gene points in the inset plot.
# fdr < 1: when enough genes are included, all adjusted p-values become 1,
#   so the gene with the minimal FDR is a random selection from all genes.
# for marking individual genes in the inset
gd <- d[fdr < 1][order(fdr)][!duplicated(gene)]

fwrite(gd, file=minpoints.csv)

# Main plot
pdf(file=main.pdf, width=3, height=3, pointsize=5)
plot(d$fdr, log='y', type='l', bty='n',
    ylab='Corrected P-value of most significant gene', ylim=c(0.01,1),
    xaxt='n', xlab='Number of genes considered\nordered from largest to smallest')
axis(side=1, at=n.genes - pretty(1:n.genes) + 1, labels=pretty(1:n.genes))
dev.off()

# Inset
pdf(file=inset.pdf, width=2, height=2, pointsize=5)
par(mar=c(2,2,0.5,0.5))
plot(d$fdr, bty='n', xlim=c(n.genes-250,n.genes+10), log='y', type='l', xaxt='n', ylim=c(0.01,1))

# we want an axis that decreases right-to-left, but R doesn't seem to support
# this natively.
xlocs <- pretty(head(1:n.genes,250))   # head() only relevant if n.genes<250
#axis(side=1, at=pretty(tail(1:n.genes,250)), labels=rev(xlocs))
axis(side=1, at=n.genes-xlocs, labels=rev(xlocs))
points(gd[,.(n.genes-top.n+1,fdr)], pch=20)
text(x=n.genes-gd$top.n+1, y=gd$fdr, labels=gd$gene, pos=1)
dev.off()

if ('snakemake' %in% ls()) {
    sink()
}
