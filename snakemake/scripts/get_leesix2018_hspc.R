#!/usr/bin/env Rscript

# detect script being run by snakemake
# if so, make a mock commandArgs function
if ('snakemake' %in% ls()) {
    logfile <- snakemake@log[[1]]
    con <- file(logfile, 'w')
    sink(con, type='output')
    sink(con, type='message')

    commandArgs <- function(...) unlist(c(
        snakemake@input['txt'],
        snakemake@output['csv'],
        snakemake@output['pdf'],
        snakemake@output['svg']
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 4) {
    stop('usage: all_mutations.txt (mendeley data download) out.csv out.pdf out.svg')
}

in.txt <- args[1]
out.csv <- args[2]
out.pdf <- args[3]
out.svg <- args[4]

for (f in c(out.csv, out.pdf, out.svg))
    if (file.exists(f))
        stop(paste0('output file ', f, ' already exists, please delete it first'))

suppressMessages(library(scan2))
suppressMessages(library(extrafont))
suppressMessages(library(svglite))
if (!("Arial" %in% fonts()))
    stop("Arial font not detected; did you load extrafonts and run font_import() with the appropriate path?")


# Read and get rid of indels
muts <- fread(in.txt)
colnames(muts)[1:4] <- c('chr', 'pos', 'refnt', 'altnt')
muts <- muts[nchar(refnt) == 1 & nchar(altnt) == 1, .(chr, pos, refnt, altnt)]
muts[, mutsig := get.3mer(muts)]

df <- as.data.frame(table(muts$mutsig))
colnames(df) <- c('MutType', 'HSPC_spectrum')
write.csv(df, file=out.csv, quote=FALSE, row.names=FALSE)

devs=list(pdf, svglite)
outs=c(out.pdf, out.svg)
for (i in 1:2) {
    devs[[i]](width=2, height=0.75, pointsize=5, file=outs[i])
    par(mar=c(0.2, 2, 2, 0.2))
    scan2::plot.sbs96(x=0, spectrum=df$HSPC_spectrum, main='HSPC spectrum (Lee-Six et al. 2018)')
    dev.off()
}


if ('snakemake' %in% ls()) {
    sink()
}
