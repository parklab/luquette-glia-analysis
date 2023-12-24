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

sens.table[, bedtype := ifelse(grepl(pattern='QBED_VERSION=', metaline), 'qbed', 'bed')]

# Each row in sens.table contains a metadata field that contains information
# about which signal is being measured. It is a semi-colon separated string with
# an arbitrary number of key=value pairs.
#
# Every metadata field contains a 'datasource' key.
sens.table[, datasource := sapply(strsplit(metaline, ';'), function(x) strsplit(x[grep('^datasource=', x)], '=')[[1]][2])]


print(table(sens.table[,.(bedtype, datasource)]))

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
    } else if (datasrc == 'cancer_snvdens') {
        key <- 'tumor'
    } else {
        stop(paste('unrecognized datasource:', datasrc))
    }
    key
}

# for non-QBEDs, there is a dataclass=X key-value pair in the metaline
sens.table[, dataclass.key := 'dataclass'] 
# For QBEDs, some manual work to figure out what signal is encoded
sens.table[bedtype == 'qbed', dataclass.key := sapply(datasource, source.to.class)]

# For QBEDs, the number of quantiles into which we divide the signal is
# stored in the metaline as QUANTILES=integer.
sens.table[bedtype == 'qbed', n.quantiles := sapply(strsplit(metaline, ';'), function(x) strsplit(x[grep('^QUANTILES=', x)], '=')[[1]][2])]
# Similarly, again for QBEDs only, the size of the tiles in the genomic
# tiling used to average the signal is stored as BINSIZE=integer.
sens.table[bedtype == 'qbed', binsize := sapply(strsplit(metaline, ';'), function(x) strsplit(x[grep('^BINSIZE=', x)], '=')[[1]][2])]

sens.table[, dataclass := mapply(function(metaline, key) strsplit(metaline[grep(paste0('^', key, '='), metaline)], '=')[[1]][2], metaline=strsplit(metaline, ';'), key=dataclass.key)]

# classify ENCODE ChIP-seq marks as either active or inactive:
#       inactive)
#           H3K27me3, H3K9me3
#       active)
#           H3K4me1, H3K4me3, H3K27ac, H3K36me3, H3K9ac
sens.table[bedtype == 'qbed' & datasource == 'encode', datasource := ifelse(dataclass == 'H3K27me3' | dataclass == 'H3K9me3', 'inactive_histone_mark', 'active_histone_mark')]

# Some manual melting:
# Each row of sens.table contains information on passtype=(pass,rescue)
# and muttype=(snv,indel) calls.  Want to split each row into 4 rows such
# that one row corresponds to exactly one [passtype,muttype] combination.
#
# static.columns DO NOT change for different [passtype,muttype] combos. Yes, that
# even includes "mean.snv.n.training*", which is the number of het germline
# SNPs used to train the AB model. Indels use the same AB model as SNVs,
# which is exclusively trained on SNP data.
static.columns <- sens.table[, .(metaline, feature, sample, width, mean.abs.gp.mu, mean.gp.sd,
                   mean.sc.dp, mean.bulk.dp, bedtype, n.quantiles, binsize,
                   datasource, dataclass.key, dataclass,
                   mean.snv.n.training, mean.snv.n.training.neighborhood)]
sens.table2 <- rbind(
    cbind(static.columns, sens.table[, .(muttype='snv', passtype='pass',
                   sens=snv.sens, n.calls=sum.snv.n.pass,
                   sum.bases.gt.sc.min.dp=sum.bases.gt.snv.sc.min.dp,
                   sum.bases.gt.bulk.min.dp=sum.bases.gt.snv.bulk.min.dp)]),
    # sens= calculation:
    # Estimate an upper bound on rescue sensitivity by assuming approximately
    # equal rates of FPs among rescued calls as among VAF-based calls.
    # For small regions where, e.g., there are no called mutations, preserve the
    # non-rescue estimate. Maybe a better idea would be to use the genome-wide
    # average rescue sensitivity.
    # Cap the maximum sensitivity at 1 since the product could exceed 1.
    cbind(static.columns, sens.table[, .(muttype='snv', passtype='rescue',
                   sens=pmin(1, snv.sens * ifelse(sum.snv.n.pass > 0 , (sum.snv.n.pass+sum.snv.n.rescue)/sum.snv.n.pass, 1)),
                   n.calls=sum.snv.n.pass + sum.snv.n.rescue,
                   sum.bases.gt.sc.min.dp=sum.bases.gt.snv.sc.min.dp,
                   sum.bases.gt.bulk.min.dp=sum.bases.gt.snv.bulk.min.dp)]),
    cbind(static.columns, sens.table[, .(muttype='indel', passtype='pass',
                   sens=indel.sens, n.calls=sum.indel.n.pass,
                   sum.bases.gt.sc.min.dp=sum.bases.gt.indel.sc.min.dp,
                   sum.bases.gt.bulk.min.dp=sum.bases.gt.indel.bulk.min.dp)]),
    cbind(static.columns, sens.table[, .(muttype='indel', passtype='rescue',
                   sens=pmin(1, indel.sens * ifelse(sum.indel.n.pass > 0, (sum.indel.n.pass+sum.indel.n.rescue)/sum.indel.n.pass, 1)),
                   n.calls=sum.indel.n.pass + sum.indel.n.rescue,
                   sum.bases.gt.sc.min.dp=sum.bases.gt.indel.sc.min.dp,
                   sum.bases.gt.bulk.min.dp=sum.bases.gt.indel.bulk.min.dp)])
)


# Get catalog-wide rates (e.g., %mutations from each sample, separate for snvs and indels)
# This is computed only on the subset of samples specified here.
metrics <- metrics[sample %in% samples]
metrics[, weight := calls.vafonly/sum(calls.vafonly), by=.(muttype)]

weights <- rbind(
    metrics[, .(passtype='pass', sample, weight=calls.vafonly/sum(calls.vafonly)), by=.(muttype)],
    metrics[, .(passtype='rescue', sample, weight=(calls.rescue+calls.vafonly)/sum(calls.rescue+calls.vafonly)), by=.(muttype)]
)

summary.table <- sens.table2[weights,,on=.(muttype, passtype, sample)][,
    .(bedtype=bedtype[1],
          # these are already in sens.table2
          #n.quantiles=n.quantiles[1],
          #binsize=binsize[1],
          group=group.tag,
          n.samples=length(samples),
          width=width[1],   # not sum(width); aggregating over samples here
          n.calls=sum(n.calls),
          weighted.sens=sum(sens*weight),
          weighted.mean.sc.dp=sum(mean.sc.dp*weight),
          weighted.mean.bulk.dp=sum(mean.bulk.dp*weight),
          weighted.mean.abs.gp.mu=sum(mean.abs.gp.mu*weight),
          weighted.mean.gp.sd=sum(mean.gp.sd*weight),
          # "snv" is always correct, even for muttype=indel
          weighted.mean.snv.n.training.neighborhood=sum(mean.snv.n.training.neighborhood*weight),
          weighted.sum.bases.gt.sc.min.dp=sum(sum.bases.gt.sc.min.dp*weight),
          weighted.sum.bases.gt.bulk.min.dp=sum(sum.bases.gt.bulk.min.dp*weight),
          weighted.fraction.bases.gt.sc.min.dp=sum(sum.bases.gt.sc.min.dp*weight)/width[1],
          weighted.fraction.bases.gt.bulk.min.dp=sum(sum.bases.gt.bulk.min.dp*weight)/width[1],
          # Weighting creates different values for snv vs. indel and pass vs. rescue
          # because the weights reflect the numbers of snv/indel/pass/rescue calls per
          # sample.  Weighted values are the mathematically accurate values to use for
          # correction, however it can cause confusion when describing general properties
          # of the genomic regions because it is natural to expect that, e.g., sequencing
          # depth of a region does not change between SNV and indel analyses.
          unweighted.sens=mean(sens),
          unweighted.mean.sc.dp=mean(mean.sc.dp),
          unweighted.mean.bulk.dp=mean(mean.bulk.dp),
          unweighted.mean.abs.gp.mu=mean(mean.abs.gp.mu),
          unweighted.mean.gp.sd=mean(mean.gp.sd),
          unweighted.mean.snv.n.training.neighborhood=mean(mean.snv.n.training.neighborhood),
          unweighted.sum.bases.gt.sc.min.dp=mean(as.numeric(sum.bases.gt.sc.min.dp)),
          unweighted.sum.bases.gt.bulk.min.dp=mean(as.numeric(sum.bases.gt.bulk.min.dp)),
          unweighted.fraction.bases.gt.sc.min.dp=mean(sum.bases.gt.sc.min.dp)/width[1],
          unweighted.fraction.bases.gt.bulk.min.dp=mean(sum.bases.gt.bulk.min.dp)/width[1]),
    by=.(datasource, dataclass, n.quantiles, binsize, feature, muttype, passtype)]

fwrite(summary.table, file=out.csv)
print(gc())

if ('snakemake' %in% ls()) {
    sink()
}
