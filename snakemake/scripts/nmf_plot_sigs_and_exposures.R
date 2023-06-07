#!/usr/bin/env Rscript

# detect script being run by snakemake
# if so, make a mock commandArgs function
if ('snakemake' %in% ls()) {
    logfile <- snakemake@log[[1]]
    con <- file(logfile, 'w')
    sink(con, type='output')
    sink(con, type='message')

    commandArgs <- function(...) unlist(c(
        snakemake@input['mat'],
        snakemake@params['muttype'],
        snakemake@output['pdf'],
        snakemake@output['expo_csv'],
        snakemake@output['proc_csv']
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 5) {
    stop("usage: plot_indel.R input.mat {SBS|ID} output.pdf out.exposures.csv out.processes.csv")
}

mat.file <- args[1]
muttype <- args[2]
outpdf <- args[3]
out.expo.csv <- args[4]
out.proc.csv <- args[5]

for (f in c(outpdf, out.expo.csv, out.proc.csv)) {
    if (file.exists(f))
        stop(paste('output file', f, 'already exists, please delete it first'))
}

suppressMessages(library(R.matlab))
suppressMessages(library(scan2))

mat <- readMat(mat.file)
names(mat$input) <- attr(mat$input, 'dimnames')[[1]]

# Get the (exposure per signature x sample) table
sample.names <- sapply(strsplit(unlist(mat$input[['sampleNames']]), split='___'), tail, 1)
expos <- t(mat$exposures)
colnames(expos) <- paste0('Signature_', 1:ncol(expos))
expo.table <- data.table(sample=sample.names, expos)

# Get the set of mutational processes (signatures) and normalize
# the order of the channels for plotting.
ctxs <- unlist(mat$input[['subtypes']])
muts <- unlist(mat$input[['types']])
if (muttype == 'SBS') {
    channel.names <- paste(ctxs, muts, sep=":")

    factorfn <- sbs96
    sigplotfn <- plot.sbs96
} else if (muttype == 'ID') {
    # Map SigProfiler ID83 channel names (like DEL_repeats_2_0) to SigProfilerMatrixGenerator
    # names (like 2:Del:R:0).
    d <- do.call(rbind, strsplit(ctxs, split="_", fixed=TRUE))
    d[,1] <- sub(pattern='DEL', replacement='Del', sub(pattern='INS', replacement="Ins", d[,1]))
    d[,2] <- sub(pattern='repeats', replacement='R', sub(pattern='MH', replacement="M", d[,2]))
    d[,3] <- sub(pattern='5+', replacement='5', d[,3], fixed=TRUE)
    d[,4] <- sub(pattern='5+', replacement='5', d[,4], fixed=TRUE)
    channel.names <- paste(d[,3], d[,1], d[,2], d[,4], sep=":")

    factorfn <- id83
    sigplotfn <- plot.id83
}

proc.table <- cbind(data.table(MutType=names(table(factorfn(c())))),
    sapply(1:ncol(mat$processes), function(i)
            setNames(mat$processes[,i], channel.names)[order(factorfn(channel.names))])
)

fwrite(expo.table, out.expo.csv)
fwrite(proc.table, out.proc.csv)

pdf(outpdf, width=18, height=2.5*ncol(mat$processes))
layout(matrix(1:(2*ncol(mat$processes)), ncol=2, byrow=TRUE), widths=c(1,3.5))
for (i in 1:ncol(mat$processes)) {
    par(mar=c(7,4,1,1))
    sigplotfn(x=0, spectrum=setNames(proc.table[[i+1]], proc.table$MutType))
    legend('topright', legend=paste0("Stability: ", round(100*mat$processStabAvg[1,i],2)), bty='n')
    par(mar=c(7,4,1,1))
    barplot(sort(setNames(expo.table[[i+1]], expo.table$sample)), las=3, border=NA, cex.names=0.80)
}

dev.off()

if ('snakemake' %in% ls()) {
    sink()
}
