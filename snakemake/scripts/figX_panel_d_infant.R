#!/usr/bin/env Rscript

# detect script being run by snakemake
# if so, make a mock commandArgs function
if ('snakemake' %in% ls()) {
    logfile <- snakemake@log[[1]]
    con <- file(logfile, 'w')
    sink(con, type='output')
    sink(con, type='message')

    commandArgs <- function(...) unlist(c(
        snakemake@input['muts'],
        snakemake@input['meta'],
        snakemake@input['cosmic'],
        snakemake@output['pdf'],
        snakemake@output['svg'],
        snakemake@output['spectra_csv'],
        snakemake@output['expo_csv']
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 7) {
    cat('muts.csv should be VAF-BASED! no rescue!\n')
    stop('usage: fig2_panel_d_infants.R muts.csv metadata.csv cosmic.csv out.pdf out.svg out.spectra.csv out.exposures.csv')
}


muts.csv <- args[1]
meta.csv <- args[2]
cosmic.csv <- args[3]
out.pdf <- args[4]
out.svg <- args[5]
out.spectra.csv <- args[6]
out.expo.csv <- args[7]

for (f in c(out.pdf, out.svg, out.spectra.csv, out.expo.csv)) {
    if (file.exists(f))
        stop(paste('output file', f, 'already exists, please delete it first'))
}

suppressMessages(library(scan2))
suppressMessages(library(pracma))
suppressMessages(library(svglite))


cosmic <- fread(cosmic.csv)
# The PTA artifact signature will be especially relevant in infants; they
# should have the highest false discovery rate.
data(snv.artifact.signature.v3)
cosmic$PTA <- snv.artifact.signature.v3

meta <- fread(meta.csv)
# These 3 subjects are the infant brains (age: 0.4, 0.6, 2.0)
muts <- fread(muts.csv)[subject %in% c(1278, 5817, 5871)]

nmuts <- muts[sample %in% meta[amp=='PTA' & type=='neuron']$sample]
omuts <- muts[sample %in% meta[amp=='PTA' & type=='oligo']$sample]
#nmuts <- fread(neuron.csv)[subject %in% c(1278, 5817, 5871)]
#omuts <- fread(oligo.csv)[subject %in% c(1278, 5817, 5871)]
# for a manual check of teenage spectra
#nmuts <- fread(neuron.csv)[subject %in% c(1465, 4638, 5559)]
#omuts <- fread(oligo.csv)[subject %in% c(1465, 4638, 5559)]

raw.neuron.spectrum <- setNames(as.vector(table(sbs96(nmuts[muttype=='snv']$mutsig))), levels(sbs96(c())))
neuron.exposures <- setNames(lsqnonneg(as.matrix(cosmic[,-1]), raw.neuron.spectrum)$x, names(cosmic)[-1])
neuron.exposures.without.artifact <- neuron.exposures
neuron.exposures.without.artifact['PTA'] <- 0
neuron.spectrum.without.artifact <- as.matrix(cosmic[,-1]) %*% neuron.exposures.without.artifact
neuron.exposures.without.artifact <- 100*neuron.exposures.without.artifact / sum(neuron.exposures.without.artifact)

raw.oligo.spectrum <- setNames(as.vector(table(sbs96(omuts[muttype=='snv']$mutsig))), levels(sbs96(c())))
oligo.exposures <- setNames(lsqnonneg(as.matrix(cosmic[,-1]), raw.oligo.spectrum)$x, names(cosmic)[-1])
oligo.exposures.without.artifact <- oligo.exposures
oligo.exposures.without.artifact['PTA'] <- 0
oligo.spectrum.without.artifact <- as.matrix(cosmic[,-1]) %*% oligo.exposures.without.artifact
oligo.exposures.without.artifact <- 100*oligo.exposures.without.artifact / sum(oligo.exposures.without.artifact)

write.csv(
    cbind(`Neuron exposure`=neuron.exposures,
          `Neuron exposure (%, PTA artifact sig removed)`=neuron.exposures.without.artifact,
          `Oligo exposure`=oligo.exposures,
          `Oligo exposure (%, PTA artifact sig removed)`=oligo.exposures.without.artifact),
    file=out.expo.csv, row.names=TRUE, quote=FALSE)

write.csv(cbind(raw.neuron.spectrum, neuron.spectrum.without.artifact,
                raw.oligo.spectrum, oligo.spectrum.without.artifact),
    file=out.spectra.csv, row.names=TRUE, quote=FALSE)

devs=list(pdf, svglite)
outs=c(out.pdf, out.svg)
for (i in 1:2) {
    devs[[i]](width=6, height=3, pointsize=5, file=outs[i])
    layout(matrix(1:6, ncol=2, byrow=TRUE))
    par(mar=c(1,4,2,1))
    plot.sbs96(x=0, spectrum=raw.neuron.spectrum, main='Infant neurons')
    plot.sbs96(x=0, spectrum=raw.oligo.spectrum, main='Infant oligos')
    plot.sbs96(x=0, spectrum=as.vector(neuron.spectrum.without.artifact), main='Infant neurons (no PTA artifact)')
    plot.sbs96(x=0, spectrum=as.vector(oligo.spectrum.without.artifact), main='Infant oligos (no PTA artifact)')
    plot.sbs96(x=0, spectrum=snv.artifact.signature.v3, main='PTA artifact signature')
    par(mar=c(4,4,0,1))
    barplot(rbind(neuron.exposures.without.artifact, oligo.exposures.without.artifact),
        beside=TRUE, las=3, col=1:2, border=NA, ylab='Percent of sSNVs from sig.')
    dev.off()
}

if ('snakemake' %in% ls()) {
    sink()
}
