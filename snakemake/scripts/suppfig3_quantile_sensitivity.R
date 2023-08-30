#!/usr/bin/env Rscript

# detect script being run by snakemake
# if so, make a mock commandArgs function
if ('snakemake' %in% ls()) {
    logfile <- snakemake@log[[1]]
    con <- file(logfile, 'w')
    sink(con, type='output')
    sink(con, type='message')

    commandArgs <- function(...) unlist(c(
        snakemake@output['csv'],
        snakemake@output['pdf'],
        snakemake@output['svg'],
        snakemake@input['sens_summary_files']
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}


args <- commandArgs(trailingOnly=T)

if (length(args) < 4)
    stop("usage: suppfig3_quantile_sensitivity.R out.csv out.pdf out.svg sens_summary1.txt [ sens_summary2.txt ... sens_summaryN.txt ]")

out.csv <- args[1]
out.pdf <- args[2]
out.svg <- args[3]
sens.files <- args[-(1:3)]

for (f in c(out.csv, out.pdf, out.svg)) {
    if (file.exists(f))
        stop(sprintf("output file %s already exists, please delete it first", f))
}

suppressMessages(library(scan2))
suppressMessages(library(svglite))

sens <- rbindlist(lapply(sens.files, fread))


n.groups <- length(unique(sens$group))
if (n.groups != length(sens.files)) {
    stop(paste('number of unique groups in table', n.groups, 'does not match number of sensitivity summary files', length(sens.files)))
}

datasources.to.plot <- c('scrnaseq', 'scatacseq', 'repliseq', 'active_histone_mark', 'inactive_histone_mark')

layout(matrix(1:(2*length(datasources.to.plot)), nrow=2, byrow=T))

for (mt in c('snv', 'indel')) {
    this.sens <- sens[muttype == mt & datasource %in% datasources.to.plot]
    # For shared ylim across plots
    ylim <- extendrange(pretty(this.sens$weighted.sens), f=0.2)
    for (dsrc in datasources.to.plot) {
        m <- dcast(this.sens[datasource==dsrc], quantile ~ dataclass+group, value.var='weighted.sens')

        # First column of m is the quantiles (e.g., 1..10)
        n.classes <- (ncol(m) - 1) /n.groups

        # try by enrichment
        mm <- m[,-1]
        #mm <- t(t(mm)/colMeans(mm))
        #ylim=c(0.8,1.2)

        par(mar=c(4,4,2,0.1))
        matplot(x=m[,1], mm, col=c('purple', 'orange', 'black','red'), lty='solid',
            lwd=1, pch=rep(letters[1:n.classes], each=n.groups), type='b', ylim=ylim,
            bty='n', xaxt=ifelse(mt == 'snv', 'n', 's'), ylab='SCAN2 sensitivity',
            xlab=ifelse(mt!='snv', 'Decile', ''), main=ifelse(mt=='snv',dsrc,''))
    }
}

if ('snakemake' %in% ls()) {
    sink()
}
