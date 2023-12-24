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
        snakemake@params['group_tag'],
        snakemake@output['csv'],
        snakemake@input['sens_files']
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}


args <- commandArgs(trailingOnly=T)

if (length(args) < 4)
    stop("usage: summarize_qbed_sensitivity.R metrics.csv group_tag output.csv sens1.txt [ sens2.txt ... sensN.txt ]")

metrics.file <- args[1]
group.tag <- args[2]
out.csv <- args[3]
sens.files <- args[-(1:3)]

for (f in c(out.csv)) {
    if (file.exists(f))
        stop(sprintf("output file %s already exists, please delete it first", f))
}

suppressMessages(library(scan2))

metrics <- fread(metrics.file)

sens.table <- rbindlist(lapply(sens.files, fread))
samples <- unique(sens.table$sample)

# Each row in sens.table contains a metadata field that contains information
# about which signal is being measured. It is a semi-colon separated string with
# an arbitrary number of key=value pairs.
#
# Every metadata field contains a 'datasource' key.
sens.table[, datasource := sapply(strsplit(metaline, ';'), function(x) strsplit(x[grep('^datasource=', x)], '=')[[1]][2])]

# Some datasources contain subclassifications (e.g., within scRNAseq there is an
# expression track for each cell type characterized by the scRNAseq). These fields
# are not uniformly named and must be handled specially here.
#   datasource:
#       encode)     mark
#       gtex)       tissue
#       repliseq)   celltype
#       scatacseq)  celltype
#       scrnaseq)   celltype
source.to.class <- function(datasrc) {
    if (datasrc == 'encode') {
        key <- 'mark'
    } else if (datasrc == 'gtex') {
        key <- 'tissue'
    } else if (datasrc == 'repliseq' | datasrc == 'scatacseq' | datasrc == 'scrnaseq') {
        key <- 'celltype'
    } else {
        stop(paste('unrecognized datasource:', datasrc))
    }
    key
}
sens.table[, dataclass.key := sapply(datasource, source.to.class)]
sens.table[, dataclass := mapply(function(metaline, key) strsplit(metaline[grep(paste0('^', key, '='), metaline)], '=')[[1]][2], metaline=strsplit(metaline, ';'), key=dataclass.key)]

# classify ENCODE ChIP-seq marks as either active or inactive:
#       inactive)
#           H3K27me3, H3K9me3
#       active)
#           H3K4me1, H3K4me3, H3K27ac, H3K36me3, H3K9ac
sens.table[datasource == 'encode', datasource := ifelse(dataclass == 'H3K27me3' | dataclass == 'H3K9me3', 'inactive_histone_mark', 'active_histone_mark')]

# now that all relevant info has been mined from metaline, remove it
sens.table$metaline <- NULL
# Order quantiles 1->10 (with support for up to 100 quantiles)
sens.table <- sens.table[feature %in% 1:100]
sens.table <- sens.table[, quantile := as.integer(feature)]
sens.table <- sens.table[, feature := paste0("quantile", quantile)]


# Get catalog-wide rates (e.g., %mutations from each sample, separate for snvs and indels)
# This is computed only on the subset of samples specified here.
metrics <- metrics[sample %in% samples]
metrics[, weight := calls.vafonly/sum(calls.vafonly), by=.(muttype)]
weights <- dcast(metrics, sample ~ muttype, value.var='weight')
colnames(weights)[2:3] <- c('snv.weight', 'indel.weight')

# Add the global weights
sens.table <- sens.table[weights,,on=.(sample)]

out.table <- sens.table[,
    .(group=group.tag,
        n.samples=length(samples),
        width=width[1],   # not sum(width); aggregating over samples here
        snv.pass=sum(sum.snv.n.pass),
        snv.rescue=sum(sum.snv.n.rescue),
        snv.unweighted.sens.vafonly=mean(snv.sens),   # just for comparison
        snv.weighted.sens.vafonly=sum(snv.sens*snv.weight),
        snv.weighted.mean.sc.dp=sum(mean.sc.dp*snv.weight),
        snv.weighted.mean.bulk.dp=sum(mean.bulk.dp*snv.weight),
        snv.weighted.mean.abs.gp.mu=sum(mean.abs.gp.mu*snv.weight),
        snv.weighted.mean.gp.sd=sum(mean.gp.sd*snv.weight),
        snv.weighted.mean.snv.n.training.neighborhood=sum(mean.snv.n.training.neighborhood*snv.weight),
        snv.weighted.sum.bases.gt.snv.sc.min.dp=sum(sum.bases.gt.snv.sc.min.dp*snv.weight),
        snv.weighted.sum.bases.gt.snv.bulk.min.dp=sum(sum.bases.gt.snv.bulk.min.dp*snv.weight),
        snv.weighted.fraction.bases.gt.snv.sc.min.dp=sum(sum.bases.gt.snv.sc.min.dp*snv.weight)/width[1],
        snv.weighted.fraction.bases.gt.snv.bulk.min.dp=sum(sum.bases.gt.snv.bulk.min.dp*snv.weight)/width[1],
        #snv.weighted.sens=sum(snv.sens*(snv.n.pass+snv.n.rescue)/sum(snv.n.pass+snv.n.rescue)),
        indel.pass=sum(sum.indel.n.pass),
        indel.rescue=sum(sum.indel.n.rescue),
        indel.unweighted.sens.vafonly=mean(indel.sens),
        indel.weighted.sens.vafonly=sum(indel.sens*indel.weight)),
        indel.weighted.mean.sc.dp=sum(mean.sc.dp*indel.weight),
        indel.weighted.mean.bulk.dp=sum(mean.bulk.dp*indel.weight),
        indel.weighted.mean.abs.gp.mu=sum(mean.abs.gp.mu*indel.weight),
        indel.weighted.mean.gp.sd=sum(mean.gp.sd*indel.weight),
        indel.weighted.mean.indel.n.training.neighborhood=sum(mean.indel.n.training.neighborhood*indel.weight),
        indel.weighted.sum.bases.gt.indel.sc.min.dp=sum(sum.bases.gt.indel.sc.min.dp*indel.weight),
        indel.weighted.sum.bases.gt.indel.bulk.min.dp=sum(sum.bases.gt.indel.bulk.min.dp*indel.weight),
        indel.weighted.fraction.bases.gt.indel.sc.min.dp=sum(sum.bases.gt.indel.sc.min.dp*indel.weight)/width[1],
        indel.weighted.fraction.bases.gt.indel.bulk.min.dp=sum(sum.bases.gt.indel.bulk.min.dp*indel.weight)/width[1],
        #indel.weighted.sens=sum(indel.sens*(indel.n.pass+indel.n.rescue)/sum(indel.n.pass+indel.n.rescue))),
    by=.(datasource, dataclass, quantile)]
fwrite(out.table, file=out.csv)

if ('snakemake' %in% ls()) {
    sink()
}
