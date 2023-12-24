#!/usr/bin/env Rscript

# detect script being run by snakemake
# if so, make a mock commandArgs function
if ('snakemake' %in% ls()) {
    logfile <- snakemake@log[[1]]
    con <- file(logfile, 'w')
    sink(con, type='output')
    sink(con, type='message')

    commandArgs <- function(...) unlist(c(
        snakemake@input['metadata'],
        snakemake@input['bulk_accessible'],
        snakemake@params['fpr_snv'],
        snakemake@params['fpr_indel'],
        snakemake@output['svg'],
        snakemake@output['pdf'],
        snakemake@output['csv'],
        snakemake@input['metrics']
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}


args <- commandArgs(trailingOnly=T)

if (length(args) < 6)
    stop(sprintf("usage: suppfigX_frac_covered.R metadata.csv bulk.accessible.txt out.svg out.pdf out.csv metrics1.txt [ metrics2.txt ... metricsN.txt ]"))

meta.file <- args[1]
bulk.accessible.file <- args[2]
fpr.snv <- as.numeric(args[3])
fpr.indel <- as.numeric(args[4])
out.svg <- args[5]
out.pdf <- args[6]
out.csv <- args[7]
metrics.files <- args[-(1:7)]

for (f in c(out.svg, out.pdf, out.csv)) {
    if (file.exists(f))
        stop(sprintf("output file %s already exists, please delete it first", f))
}

# DO WORK
suppressMessages(library(data.table))
suppressMessages(library(svglite))

meta <- fread(meta.file)
meta$ageclass <- factor(meta$ageclass, levels=c('infant', 'adolescent', 'adult', 'elderly'), ordered=T)

bulk.accessible <- fread(bulk.accessible.file)
accessible.bp <- sum(bulk.accessible$n.basepairs)

samples <- list(`PTA neuron`=meta[amp=='PTA' & type=='neuron' & outlier == 'NORMAL']$sample,
    `PTA oligo`=meta[amp=='PTA' & type=='oligo' & outlier == 'NORMAL']$sample,
    `MDA oligo`=meta[amp=='MDA' & type=='oligo' & outlier == 'NORMAL']$sample,
    `MDA GFAP`=meta[amp=='MDA' & type=='mixed' & outlier == 'NORMAL']$sample)

metrics <- rbindlist(lapply(metrics.files, fread))
# pmin: the estimated number of FPs should not exceed the actual number of calls
# Note that a single FPR for all mutation burdens is used, but this only an
# approximation. A more ideal estimate would adjust FPR for mutation burden,
# since cells with fewer true mutations use less permissive calling parameters
# and thus likely commit fewer FPs.
metrics[muttype == 'snv', fps := pmin(calls.vafonly+calls.rescue, analyzable.basepairs/1e6 * fpr.snv)]
metrics[muttype == 'indel', fps := pmin(calls.vafonly+calls.rescue, analyzable.basepairs/1e6 * fpr.indel)]
metrics[, fdr := ifelse(calls.vafonly+calls.rescue == 0, 0, fps / (calls.vafonly+calls.rescue))]
fwrite(metrics, file=out.csv)

# For both mean.deth and mapd, muttype=snv or indel does not matter. The same value
# is used for each.
mean.depth <- lapply(samples, function(s) metrics[muttype=='snv' & sample %in% s]$mean.depth)
mapd <- lapply(samples, function(s) metrics[muttype=='snv' & sample %in% s]$mapd)

analyzable.snvs <- lapply(samples, function(s) metrics[muttype=='snv' & sample %in% s]$analyzable.basepairs / accessible.bp)
analyzable.indels <- lapply(samples, function(s) metrics[muttype=='indel' & sample %in% s]$analyzable.basepairs / accessible.bp)

# All outlier samples have scaling factor=NA
burden.scaling.snvs <- lapply(samples, function(s) metrics[muttype=='snv' & sample %in% s & !is.na(burden.scaling.factor)]$burden.scaling.factor)
burden.scaling.indels <- lapply(samples, function(s) metrics[muttype=='indel' & sample %in% s & !is.na(burden.scaling.factor)]$burden.scaling.factor)

vafonly.sens.snvs <- lapply(samples, function(s) metrics[muttype=='snv' & sample %in% s]$vafonly.sensitivity)
vafonly.sens.indels <- lapply(samples, function(s) metrics[muttype=='indel' & sample %in% s]$vafonly.sensitivity)


# Must remove calls=0 cases since this leads to division by 0 when
# estimating sensitivity increase by rescue.  Rescue sens is estimated
# by comparing (#calls rescued or passed) / (#calls passed)
rescue.sens.ub.snvs <- lapply(samples, function(s) metrics[muttype=='snv' & sample %in% s & calls.vafonly > 0]$rescue.sensitivity.upper.bound)
rescue.sens.ub.indels <- lapply(samples, function(s) metrics[muttype=='indel' & sample %in% s & calls.vafonly > 0]$rescue.sensitivity.upper.bound)

fdr.snvs <- lapply(samples, function(s) metrics[muttype=='snv' & sample %in% s]$fdr)
fdr.indels <- lapply(samples, function(s) metrics[muttype=='indel' & sample %in% s]$fdr)

weighted.fdr.snvs <- sapply(samples, function(s) metrics[muttype=='snv' & sample %in% s & !is.na(fps), sum(fps)/sum(calls.vafonly + calls.rescue)])
weighted.fdr.indels <- sapply(samples, function(s) metrics[muttype=='indel' & sample %in% s & !is.na(fps), sum(fps)/sum(calls.vafonly + calls.rescue)])

print(unlist(weighted.fdr.snvs))
print(unlist(weighted.fdr.indels))


box.and.strip <- function(x, ylim=range(x), ...) {
print(x)
#print(ylab)
    boxplot(x,
        las=3,
        col=c('#444444', '#DE556E', 'orange', 'purple'),
        outline=F, pars=list(whisklty='solid', bty='n'),
        ylim=ylim, ...)

    stripchart(x, vertical=TRUE, pch=20, add=TRUE, method='jitter')
}

devs <- list(svglite, pdf)
outs <- c(out.svg, out.pdf)
for (i in 1:2) {
    ncols=2
    nrows=10
    devs[[i]](file=outs[i], width=ncols*1.75, height=nrows*2.5)#, pointsize=5)
    layout(matrix(1:(ncols*nrows), ncol=ncols, byrow=TRUE))

    par(mar=c(6,4,1,1))
    box.and.strip(mean.depth, ylab='Mean sequencing depth')
    box.and.strip(mapd, ylab='MAPD')

    box.and.strip(analyzable.snvs, ylim=range(c(unlist(analyzable.snvs), unlist(analyzable.indels))),
        ylab='Fraction of genome >= min. depth')
    box.and.strip(analyzable.indels, ylim=range(c(unlist(analyzable.snvs), unlist(analyzable.indels))),
        ylab='', yaxt='n')

    box.and.strip(burden.scaling.snvs, ylim=range(c(unlist(burden.scaling.snvs), unlist(burden.scaling.indels))),
        ylab='Genome-wide burden scaling factor')
    box.and.strip(burden.scaling.indels, ylim=range(c(unlist(burden.scaling.snvs), unlist(burden.scaling.indels))),
        ylab='', yaxt='n')


    box.and.strip(vafonly.sens.snvs, ylim=range(c(unlist(vafonly.sens.snvs), unlist(vafonly.sens.indels))),
        ylab='SCAN2 sensitivity (VAF-based)')
    box.and.strip(vafonly.sens.indels, ylim=range(c(unlist(vafonly.sens.snvs), unlist(vafonly.sens.indels))),
        ylab='', yaxt='n')

    box.and.strip(rescue.sens.ub.snvs, ylim=range(c(unlist(rescue.sens.ub.snvs), unlist(rescue.sens.ub.indels))),
        ylab='SCAN2 rescue sensitivity upper bound')
    box.and.strip(rescue.sens.ub.indels, ylim=range(c(unlist(rescue.sens.ub.snvs), unlist(rescue.sens.ub.indels))),
        ylab='', yaxt='n')

    box.and.strip(fdr.snvs, ylim=range(c(unlist(fdr.snvs), unlist(fdr.indels))),
        ylab='SCAN2 false discovery rate')
    box.and.strip(fdr.indels, ylim=range(c(unlist(fdr.snvs), unlist(fdr.indels))),
        ylab='', yaxt='n')


    fdr <- meta[metrics, , on=.(sample)]
    s.fdr <- fdr[amp=='PTA' & muttype=='snv' & outlier=='NORMAL', .(age,fdr,type,color)][order(age)]
    plot(s.fdr[,.(age, fdr)], pch=20, col=s.fdr$color, ylim=0:1, ylab='False discovery rate', xlab='Age', bty='n')
    lines(s.fdr[type=='neuron']$age, predict(loess(fdr~age, data=s.fdr[type=='neuron'])), lwd=2, col='#444444')
    lines(s.fdr[type=='oligo']$age, predict(loess(fdr~age, data=s.fdr[type=='oligo'])), lwd=2, col='#DE556E')

    i.fdr <- fdr[amp=='PTA' & muttype=='indel' & outlier=='NORMAL', .(age,fdr,type,color)][order(age)]
    plot(i.fdr[,.(age, fdr)], pch=20, col=i.fdr$color, ylim=0:1, ylab='', yaxt='n', xlab='Age', bty='n')
    lines(i.fdr[type=='neuron']$age, predict(loess(fdr~age, data=i.fdr[type=='neuron'])), lwd=2, col='#444444')
    lines(i.fdr[type=='oligo']$age, predict(loess(fdr~age, data=i.fdr[type=='oligo'])), lwd=2, col='#DE556E')
    

    n.snv <- fdr[amp=='PTA' & muttype=='snv', .(n.calls=sum(calls.vafonly + calls.rescue)), by=.(ageclass,type)][order(ageclass)]
    n.snv.matrix <- rbind(n.snv[type=='neuron']$n.calls, n.snv[type=='oligo']$n.calls)
    rownames(n.snv.matrix) <- c('neuron', 'oligo')
    colnames(n.snv.matrix) <- n.snv[type=='neuron']$ageclass
print(n.snv.matrix)

    n.indel <- fdr[amp=='PTA' & muttype=='indel', .(n.calls=sum(calls.vafonly + calls.rescue)), by=.(ageclass,type)][order(ageclass)]
    n.indel.matrix <- rbind(n.indel[type=='neuron']$n.calls, n.indel[type=='oligo']$n.calls)
    rownames(n.indel.matrix) <- c('neuron', 'oligo')
    colnames(n.indel.matrix) <- n.indel[type=='neuron']$ageclass
print(n.indel.matrix)

    barplot(n.snv.matrix,
        beside=TRUE, border=NA, ylab='SNVs in catalog',
        col=c('#444444', '#DE556E'), las=3, ylim=c(0, max(pretty(n.snv.matrix))))
    barplot(n.indel.matrix,
        beside=TRUE, border=NA, ylab='Indels in catalog',
        col=c('#444444', '#DE556E'), las=3, ylim=c(0, max(pretty(n.indel.matrix))))

    bp <- barplot(100*weighted.fdr.snvs,
        las=3, ylim=c(0,10),
        col=c('#444444', '#DE556E', 'orange', 'purple'),
        ylab='Mutation catalog false positives (%)')
    text(x=bp, y=100*weighted.fdr.snvs, labels=round(100*weighted.fdr.snvs,1), pos=3)
    bp <- barplot(100*weighted.fdr.indels,
        las=3, ylim=c(0,10),
        col=c('#444444', '#DE556E', 'orange', 'purple'),
        ylab='', yaxt='n')
    text(x=bp, y=100*weighted.fdr.indels, labels=round(100*weighted.fdr.indels,1), pos=3)

    plot(1, pch=NA, xaxt='n', xlab='', yaxt='n', ylab='', bty='n')
    legend('center', bty='n',
        fill=c('#444444', '#DE556E', 'orange', 'purple'),
        legend=c('PTA neurons', 'PTA oligos', 'MDA oligos', 'MDA GFAP'))

    dev.off()
}

if ('snakemake' %in% ls()) {
    sink()
}
