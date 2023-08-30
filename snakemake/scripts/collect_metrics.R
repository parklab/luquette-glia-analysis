#!/usr/bin/env Rscript

# detect script being run by snakemake
# if so, make a mock commandArgs function
if ('snakemake' %in% ls()) {
    logfile <- snakemake@log[[1]]
    con <- file(logfile, 'w')
    sink(con, type='output')
    sink(con, type='message')

    commandArgs <- function(...) unlist(c(
        snakemake@input['full_object'],
        snakemake@input['mapd'],
        snakemake@output['metrics'],
        snakemake@output['dp_cdf']
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}


args <- commandArgs(trailingOnly=T)

if (length(args) != 4)
    stop(sprintf("usage: collect_metrics.R scan2_full_object.rda mapd.txt output_metrics.csv depth_cdf.csv"))

in.rda <- args[1]
in.mapd <- args[2]
out.metrics <- args[3]
out.depth <- args[4]

for (f in c(out.metrics, out.depth)) {
    if (file.exists(f))
        stop(sprintf("output file %s already exists, please delete it first", f))
}

suppressMessages(library(scan2))

mapd <- scan(in.mapd, what=0)

load(in.rda, verb=TRUE)

dp.cdf <- data.table(dp=0:(nrow(results@depth.profile$dptab)-1), n.basepairs=rowSums(results@depth.profile$dptab))

call.metrics <- function(object, mt) {
    gatk <- object@gatk

    calls.vafonly <- gatk[muttype == mt, sum(pass, na.rm=TRUE)]
    calls.rescue <- gatk[muttype == mt, sum(rescue, na.rm=TRUE)]
    vafonly.sensitivity <- gatk[muttype == mt & resampled.training.site==TRUE, mean(training.pass, na.rm=TRUE)]
    rescue.sensitivity.upper.bound.factor <- gatk[muttype == mt, sum(pass | rescue, na.rm=TRUE)/sum(pass, na.rm=TRUE)]

    # row 2 is the middle 50% of sequencing depth, which is the portion of
    # the genome we use to extrapolate total burden. this trimmed mean
    # approach, which removes anomolously low and high depth genomic regions,
    # helps to stabilize the extrapolation.
    mb <- object@mutburden[[mt]][2,]
    # dividing by calls.vafonly: we extrapolate genome-wide burden using only
    # the calls in the middle 50% (mb$ncalls), but there is a desire to have a
    # scaling factor that is applicable to the total number of called mutations
    burden.scaling.factor.per.gbp <- mb$ncalls / mb$callable.sens / mb$callable.bp * 1e9 / 2 / calls.vafonly

    # to calculate number of bases passsing min depth filters
    sfp <- object@static.filter.params[[mt]]
    dptab <- object@depth.profile$dptab
    analyzable.basepairs <- sum(dptab[(sfp$min.sc.dp+1):nrow(dptab),(sfp$min.bulk.dp+1):ncol(dptab)])

    data.table(
        analyzable.basepairs=analyzable.basepairs,
        calls.vafonly=calls.vafonly,
        calls.rescue=calls.rescue,
        vafonly.sensitivity=vafonly.sensitivity,
        rescue.sensitivity.upper.bound=rescue.sensitivity.upper.bound.factor * vafonly.sensitivity,
        burden.scaling.factor.per.gbp=burden.scaling.factor.per.gbp,
        burden.scaling.factor=burden.scaling.factor.per.gbp * get.gbp.by.genome(object)
    )
}

metrics <- rbind(
    data.table(sample=results@single.cell, mapd=mapd, mean.depth=mean.coverage(results)['single.cell'], muttype='snv', call.metrics(results, 'snv')),
    data.table(sample=results@single.cell, mapd=mapd, mean.depth=mean.coverage(results)['single.cell'], muttype='indel', call.metrics(results, 'indel'))
)

fwrite(dp.cdf, file=out.depth)
fwrite(metrics, file=out.metrics)

if ('snakemake' %in% ls()) {
    sink()
}
