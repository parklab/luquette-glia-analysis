#!/usr/bin/env Rscript

# detect script being run by snakemake
# if so, make a mock commandArgs function
if ('snakemake' %in% ls()) {
    logfile <- snakemake@log[[1]]
    con <- file(logfile, 'w')
    sink(con, type='output')
    sink(con, type='message')

    commandArgs <- function(...) {
        ret <- unlist(c(
            snakemake@output['pdf'],
            snakemake@output['svg'],
            snakemake@output['csv'],
            paste0(snakemake@params['tag1'], '=', snakemake@input['cors1'])))
        if ('tag2' %in% names(snakemake@params))
            ret <- c(ret, paste0(snakemake@params['tag2'], '=', snakemake@input['cors2']))
        ret
    }
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 4 & length(args) != 5) {
    stop('usage: fig6_panel_a.R out.pdf out.svg out.csv tag1=qbed_cors1.csv [ tag2=qbed_cors2.csv ]')
}


out.pdf <- args[1]
out.svg <- args[2]
out.csv <- args[3]
cors1.file <- args[4]

do.two <- FALSE
if (length(args) == 5) {
    do.two <- TRUE
    cors2.file <- args[5]
}

for (f in c(out.pdf, out.svg, out.csv)) {
    if (file.exists(f))
        stop(paste('output file', f, 'already exists, please delete it first'))
}

suppressMessages(library(data.table))
suppressMessages(library(svglite))

cors1.tag <- strsplit(cors1.file, '=')[[1]][1]
cors1.file <- strsplit(cors1.file, '=')[[1]][2]
print(cors1.tag)
print(cors1.file)
npv <- fread(cors1.file)[order(Correlation)]
npv$Cancer <- sapply(npv$Signal, function(signal) mutenrich::read.bed.metadata(string=signal, is.qbed=TRUE)[['tumor']])
npv$Signal <- NULL
npv[, group := cors1.tag]

all.pv <- npv
if (do.two) {
    cors2.tag <- strsplit(cors2.file, '=')[[1]][1]
    cors2.file <- strsplit(cors2.file, '=')[[1]][2]
    print(cors2.tag)
    print(cors2.file)
    opv <- fread(cors2.file)[order(Correlation)]
    opv$Cancer <- sapply(opv$Signal, function(signal) mutenrich::read.bed.metadata(string=signal, is.qbed=TRUE)[['tumor']])
    opv$Signal <- NULL
    opv[, group := cors2.tag]

    all.pv <- rbind(all.pv, opv)
}

fwrite(all.pv, file=out.csv)

colors <- setNames(rep('grey',37), npv$Cancer)
colors[c('CNS-GBM', 'CNS-Medullo', 'CNS-Oligo', 'CNS-PiloAstro')] <-
    c('orange','black','red','purple')

devs=list(pdf, svglite)
outs=c(out.pdf, out.svg)
for (i in 1:2) {
    if (do.two) {
        devs[[i]](width=6.5, height=2, pointsize=5, file=outs[i])
        layout(t(1:2))
    } else {
        # decrease width
        devs[[i]](width=3.25, height=2, pointsize=5, file=outs[i])
    }
        
    par(mar=c(8,4,3,1))
    if (do.two) {
        barplot(opv$Correlation, col=colors[opv$Cancer], border=F, las=3, ylim=c(-0.03,0.40),
            cex.names=0.8, names.arg=opv$Cancer, ylab='Correlation')
    }
    barplot(npv$Correlation, col=colors[npv$Cancer], border=F, las=3, ylim=c(-0.03,0.40),
        cex.names=0.8, names.arg=npv$Cancer, ylab='Correlation')
    dev.off()
}

if ('snakemake' %in% ls()) {
    sink()
}
