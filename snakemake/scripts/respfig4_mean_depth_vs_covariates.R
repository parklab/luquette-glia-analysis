#!/usr/bin/env Rscript

# detect script being run by snakemake
# if so, make a mock commandArgs function
if ('snakemake' %in% ls()) {
    logfile <- snakemake@log[[1]]
    con <- file(logfile, 'w')
    sink(con, type='output')
    sink(con, type='message')

    commandArgs <- function(...) unlist(c(
        snakemake@input['metrics'],
        snakemake@input['meta'],
        snakemake@output['pdf'],
        snakemake@output['csv']
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}


args <- commandArgs(trailingOnly=T)

if (length(args) != 4)
    stop("usage: respfig4_mean_depth_vs_coviarates.R in.metrics.csv in.metadata.csv out.pdf out.csv")

in.metrics <- args[1]
in.meta <- args[2]
out.pdf <- args[3]
out.csv <- args[4]

for (f in c(out.pdf, out.csv)) {
    if (file.exists(f))
        stop(sprintf("output file %s already exists, please delete it first", f))
}

# DO WORK
library(data.table)

x <- fread(in.metrics)
meta <- fread(in.meta)

x <- meta[x,,on=.(sample)]
x <- x[amp=='PTA'] # & mean.depth<48]

x[, burden := calls.vafonly*burden.scaling.factor]
x[, calls.all := calls.vafonly+calls.rescue]

fwrite(x[,.(sample,donor,age,amp,selection,outlier,sex,muttype,mean.depth,mapd,burden.scaling.factor,calls.all,burden)], file=out.csv)
pdf(file=out.pdf, pointsize=7, width=7, height=3.5)
layout(matrix(1:8, nrow=2))
for (col in c('mapd', 'burden.scaling.factor', 'calls.all', 'burden')) {
    par(mar=c(2,2,2,1))
    xx <- x[muttype=='snv']
    plot(xx$mean.depth, xx[[col]], col=xx$color, pch=ifelse(xx$outlier == "NORMAL", 20, 4), bty='l', xlab='', main=paste('snv', col), ylab='')
    if (col == 'burden') {
        legend('topright', legend=c('Neuron', 'Oligo', 'Outlier'), pch=c(20,20,4),
            col=c(x[type=='neuron']$color[1], x[type=='oligo']$color[2], 'darkgrey'), cex=1.5, bty='n')
    }
    abline(v=47, col='grey', lty='dotted')
    xx <- x[muttype=='indel']
    plot(xx$mean.depth, xx[[col]], col=xx$color, pch=ifelse(xx$outlier == "NORMAL", 20, 4), bty='l', xlab='', main=paste('indel', col), ylab='')
    abline(v=47, col='grey', lty='dotted')
}
dev.off()

if ('snakemake' %in% ls()) {
    sink()
}
