#!/usr/bin/env Rscript

# detect script being run by snakemake
# if so, make a mock commandArgs function
if ('snakemake' %in% ls()) {
    logfile <- snakemake@log[[1]]
    con <- file(logfile, 'w')
    sink(con, type='output')
    sink(con, type='message')

    commandArgs <- function(...) unlist(c(
        snakemake@params['group_tag'],
        snakemake@output['csv'],
        snakemake@input['sens_files']
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}


args <- commandArgs(trailingOnly=T)

if (length(args) < 3) {
    cat("This script only produces weighted sensitivity for VAF-based SCAN2 calls\n")
    stop("usage: summarize_qbed_sensitivity.R group_tag output.csv sens1.txt [ sens2.txt ... sensN.txt ]")
}

group.tag <- args[1]
out.csv <- args[2]
sens.files <- args[-(1:2)]

for (f in c(out.csv)) {
    if (file.exists(f))
        stop(sprintf("output file %s already exists, please delete it first", f))
}

suppressMessages(library(scan2))

sens.tables <- lapply(sens.files, fread)
samples <- sapply(sens.tables, function(t) unique(t$sample))  # should only be 1 sample value per table
cat('got', length(samples), 'samples:\n')
print(samples)

# sanity check: all tables are parallel tilings
first.sens <- sens.tables[[1]]
for (i in 2:length(sens.tables)) {
    t <- sens.tables[[i]]
    if (any(first.sens$chr != t$chr))
        stop(paste('column CHR in file', i, '=', sens.files[i], 'does not match first file. all files must match'))
    if (any(first.sens$start != t$start))
        stop(paste('column START in file', i, '=', sens.files[i], 'does not match first file. all files must match'))
    if (any(first.sens$end != t$end))
        stop(paste('column END in file', i, '=', sens.files[i], 'does not match first file. all files must match'))
}

snv.counts <- setNames(sapply(sens.tables, function(t) sum(t$sum.snv.pass.n.calls)), samples)
snv.weights <- snv.counts / sum(snv.counts)
indel.counts <- setNames(sapply(sens.tables, function(t) sum(t$sum.indel.pass.n.calls)), samples)
indel.weights <- indel.counts / sum(indel.counts)

print(snv.counts)
print(indel.counts)

final.table <- first.sens[,.(chr,start,end,width,
    group.tag=group.tag, n.samples=length(samples),
    total.snvs=sum(snv.counts), total.indels=sum(indel.counts))]

str(sapply(sens.tables, function(t) t$snv.sens))
str(colMeans(sapply(sens.tables, function(t) t$snv.sens)))

final.table$unweighted.snv.sens <- rowMeans(sapply(sens.tables, function(t) t$snv.sens))
final.table$weighted.snv.sens <- rowSums(sapply(sens.tables, function(t) t$snv.sens*snv.weights[t$sample[1]]))
final.table$unweighted.indel.sens <- rowMeans(sapply(sens.tables, function(t) t$indel.sens))
final.table$weighted.indel.sens <- rowSums(sapply(sens.tables, function(t) t$indel.sens*indel.weights[t$sample[1]]))

fwrite(final.table, file=out.csv)
print(gc())

if ('snakemake' %in% ls()) {
    sink()
}
