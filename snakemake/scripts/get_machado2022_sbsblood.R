#!/usr/bin/env Rscript

# detect script being run by snakemake
# if so, make a mock commandArgs function
if ('snakemake' %in% ls()) {
    logfile <- snakemake@log[[1]]
    con <- file(logfile, 'w')
    sink(con, type='output')
    sink(con, type='message')

    commandArgs <- function(...) unlist(c(
        snakemake@input['xlsx'],
        snakemake@params['sheet_number'],
        snakemake@output['supptab8'],
        snakemake@output['csv'],
        snakemake@output['pdf'],
        snakemake@output['svg']
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 6) {
    stop('usage: all_supplementary_tables.xlsx sheet_number (=8) out.supptab8.csv out.sbsblood.csv out.pdf out.svg')
}

in.xlsx <- args[1]
sheet.number <- as.integer(args[2])
out.supptab8.csv <- args[3]
out.sbsblood.csv <- args[4]
out.pdf <- args[5]
out.svg <- args[6]

for (f in c(out.sbsblood.csv, out.pdf, out.svg))
    if (file.exists(f))
        stop(paste0('output file ', f, ' already exists, please delete it first'))

suppressMessages(library(scan2))
suppressMessages(library(pracma))
suppressMessages(library(extrafont))
suppressMessages(library(svglite))
suppressMessages(library(rio))
if (!("Arial" %in% fonts()))
    stop("Arial font not detected; did you load extrafonts and run font_import() with the appropriate path?")


rio::convert(in_file=in.xlsx, out_file=out.supptab8.csv, in_opts=list(sheet=sheet.number))
tmp <- read.csv(out.supptab8.csv)
tmp <- tmp[tmp[,1] %in% c('Signature', 'SBSblood'),]
df <- data.frame(
    MutType=paste0(substr(tmp[1,-1], 1, 3), ':', substr(tmp[1,-1],2,2), '>', substr(tmp[1,-1],6,6)),
    SBSblood=as.numeric(tmp[2,-1]))

write.csv(df, file=out.sbsblood.csv, quote=FALSE, row.names=FALSE)

devs=list(pdf, svglite)
outs=c(out.pdf, out.svg)
for (i in 1:2) {
    devs[[i]](width=2, height=0.75, pointsize=5, file=outs[i])
    par(mar=c(0.2, 2, 2, 0.2))
    scan2::plot.sbs96(x=0, spectrum=df$SBSblood, main='SBSblood (Machado et al. 2022)')
    dev.off()
}


if ('snakemake' %in% ls()) {
    sink()
}
