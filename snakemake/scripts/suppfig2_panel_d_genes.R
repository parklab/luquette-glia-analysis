#!/usr/bin/env Rscript

# detect script being run by snakemake
# if so, make a mock commandArgs function
if ('snakemake' %in% ls()) {
    logfile <- snakemake@log[[1]]
    con <- file(logfile, 'w')
    sink(con, type='output')
    sink(con, type='message')

    commandArgs <- function(...) unlist(c(
        snakemake@input['neuron'],
        snakemake@input['oligo'],
        snakemake@output['csv'],
        snakemake@output['pdf'],
        snakemake@output['svg'],
        snakemake@output['jpeg']
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 6) {
    stop('usage: fig1_genes_suppl.R neuron.SUMMARY.rda oligo.SUMMARY.rda out.csv out.pdf out.svg out.jpeg')
}


neuron.rda <- args[1]
oligo.rda <- args[2]
out.csv <- args[3]
out.pdf <- args[4]
out.svg <- args[5]
out.jpeg <- args[6]

for (f in c(out.csv, out.pdf, out.svg, out.jpeg)) {
    if (file.exists(f))
        stop(paste('output file', f, 'already exists, please delete it first'))
}

suppressMessages(library(scan2))
suppressMessages(library(svglite))


load(neuron.rda)
neuron <- es[[1]]
neuron$gene <- rownames(neuron)
setDT(neuron)
neuron <- neuron[, .(gene, pval, enr, obs, perm.mean, enr.perm.lb.95, enr.perm.ub.95)]
colnames(neuron)[-1] <- paste0('neuron.', colnames(neuron)[-1])
setkey(neuron, gene)

load(oligo.rda)
oligo <- es[[1]]
oligo$gene <- rownames(oligo)
setDT(oligo)
oligo <- oligo[, .(gene, pval, enr, obs, perm.mean, enr.perm.lb.95, enr.perm.ub.95)]
colnames(oligo)[-1] <- paste0('oligo.', colnames(oligo)[-1])
setkey(oligo, gene)

combined <- neuron[oligo][order(pmin(neuron.pval, oligo.pval))]
# outside is the special feature created by bedenrich.r to contain all
# mutations that don't intersect the feature set.
combined <- combined[gene != 'outside']

fwrite(combined, file=out.csv)

# just for plotting. don't save this in the table
minp <- pmin(combined$neuron.pval, combined$oligo.pval)
xlim <- combined[, range(pretty(range(-log10(neuron.pval) * sign(log(neuron.enr)), na.rm=TRUE)))]
ylim <- combined[, range(pretty(range(-log10(oligo.pval) * sign(log(oligo.enr)), na.rm=TRUE)))]

devs=list(pdf, svglite, jpeg)
outs=c(out.pdf, out.svg, out.jpeg)
for (i in 1:3) {
    if (i == 3)
        devs[[i]](width=4, height=4, pointsize=5, file=outs[i], unit='in', res=600)
    else
        devs[[i]](width=4, height=4, pointsize=5, file=outs[i])
    par(mar=c(4,4,3,1))

    plot(combined[,.(-log10(neuron.pval) * sign(log(neuron.enr)), -log10(oligo.pval) * sign(log(oligo.enr)))],
        pch=20, cex=1/2,
        xlab='Neuron -log10(enrichment P-value)',
        ylab='Oligo -log10(enrichment P-value)',
        main='Negative value: under-mutated, Positive value: over-mutated',
        col=ifelse(minp < 0.05, 'black', 'grey'), xlim=xlim, ylim=ylim)

    basicPlotteR::addTextLabels(xCoords=combined[minp < 0.05, -log10(neuron.pval)*sign(log(neuron.enr))],
                                yCoords=combined[minp < 0.05, -log10(oligo.pval)*sign(log(oligo.enr))],
        labels=combined[minp < 0.05]$gene, col.line=1, col.label=1, cex.label=0.9)

    dev.off()
}

if ('snakemake' %in% ls()) {
    sink()
}
