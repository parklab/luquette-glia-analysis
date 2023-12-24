#!/usr/bin/env Rscript

# detect script being run by snakemake
# if so, make a mock commandArgs function
if ('snakemake' %in% ls()) {
    logfile <- snakemake@log[[1]]
    con <- file(logfile, 'w')
    sink(con, type='output')
    sink(con, type='message')

    commandArgs <- function(...) unlist(c(
        snakemake@input['cosmic'],
        snakemake@input['sigb'],
        snakemake@output['pdf'],
        snakemake@output['suppl_pdf'],
        snakemake@output['svg'],
        snakemake@output['suppl_svg'],
        snakemake@output['spectra_csv'],
        snakemake@output['expo_pdf'],
        snakemake@output['expo_svg'],
        snakemake@output['expo_csv'],
        snakemake@output['spectrum_early_vs_late_pdf'],
        snakemake@output['spectrum_early_vs_late_csv'],
        snakemake@input[['muts']]   # list
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 12) {
    cat('muts.csv must be a combined 2-sample mutation table as produced by fit_mrca.R\n')
    cat('IMPORTANT: THIS SCRIPT IS NOT GENERALIZABLE! One of the 3 OL pairs has a large number of shared sSNVs, meaning its MRCA was later in life. The mutation spectrum of that pair is plotted against the other two pairs to contrast the mostly-developmental mutations from the combination of aging and developmental mutations.\n')
    stop('usage: figX_panels_bc.R cosmic.csv sigb.csv out.pdf out.supplement.pdf out.svg out.supplement.svg out.spectra.csv out.exposures_barplot.pdf out.exposures_barplot.svg out.exposures_barplot.csv out_early_late.pdf out_early_late.csv muts1.csv [ muts2.csv ... mutsN.csv ]')
}


cosmic.csv <- args[1]
sigb.csv <- args[2]
out.pdf <- args[3]
out.suppl.pdf <- args[4]
out.svg <- args[5]
out.suppl.svg <- args[6]
out.spectra.csv <- args[7]
out.expo.pdf <- args[8]
out.expo.svg <- args[9]
out.expo.csv <- args[10]
out.early_vs_late.pdf <- args[11]
out.early_vs_late.csv <- args[12]
mut.csvs <- args[-(1:12)]

for (f in c(out.pdf, out.suppl.pdf, out.svg, out.suppl.svg, out.spectra.csv, out.expo.pdf, out.expo.svg, out.expo.csv, out.early_vs_late.pdf, out.early_vs_late.csv)) {
    if (file.exists(f))
        stop(paste('output file', f, 'already exists, please delete it first'))
}

suppressMessages(library(scan2))
suppressMessages(library(pracma))
suppressMessages(library(svglite))


cosmic <- fread(cosmic.csv)
lod <- fread(sigb.csv)
cosmic$SigB <- lod$B
cosmic$PTA <- get(data(snv.artifact.signature.v3))   # loads `snv.artifact.signature.v3`

mut.list <- lapply(mut.csvs, function(f) fread(f)[,.(mutsig,shared,private)])

late.pair.idx <- which.max(sapply(mut.list, function(d) nrow(d[shared==TRUE])))
muts.late.pair <- mut.list[[late.pair.idx]]
muts.early.pairs <- do.call(rbind, mut.list[-late.pair.idx])

muts <- do.call(rbind, mut.list)

raw.shared.spectrum <- as.vector(table(sbs96(muts[shared==TRUE]$mutsig)))
sh <- setNames(lsqnonneg(as.matrix(cosmic[,-1]), raw.shared.spectrum)$x,
               colnames(cosmic)[-1])
sh.sbs1.pct <- sh['SBS1'] / sum(sh[names(sh) != 'SigB'])
# SigB has extremely low contribution to the shared spectrum, but do this for consistency
tmp.sh <- sh
tmp.sh['SigB'] <- 0 # exclude signature B contribution from plotted spectrum
tmp.sh['PTA'] <- 0 # exclude PTA contribution from plotted spectrum
sh.expo.pct <- tmp.sh/sum(tmp.sh)*100
reconstructed.shared.spectrum <- setNames(as.vector(as.matrix(cosmic[,-1]) %*% tmp.sh),
    as.vector(cosmic$MutType))

raw.private.spectrum <- as.vector(table(sbs96(muts[private==TRUE]$mutsig)))
pr <- setNames(lsqnonneg(as.matrix(cosmic[,-1]), raw.private.spectrum)$x,
               colnames(cosmic)[-1])
pr.sbs1.pct <- pr['SBS1'] / sum(pr[names(pr) != 'SigB'])
tmp.pr <- pr
tmp.pr['SigB'] <- 0 # exclude signature B contribution from plotted spectrum
tmp.pr['PTA'] <- 0
pr.expo.pct <- tmp.pr/sum(tmp.pr)*100
reconstructed.private.spectrum <- setNames(as.vector(as.matrix(cosmic[,-1]) %*% tmp.pr),
    as.vector(cosmic$MutType))

write.csv(
    cbind(`Shared sSNV exposure`=sh,
        `Shared sSNV exposure (%, sig B and PTA removed)`=sh.expo.pct,
        `Private sSNV exposure`=pr,
        `Private sSNV exposure (%, sig B and PTA removed)`=pr.expo.pct),
    file=out.expo.csv, row.names=TRUE, quote=FALSE)

write.csv(cbind(raw.shared.spectrum, reconstructed.shared.spectrum,
                raw.private.spectrum, reconstructed.private.spectrum),
    file=out.spectra.csv, row.names=TRUE, quote=FALSE)

# reconstructed spectra
devs=list(pdf, svglite)
outs=c(out.pdf, out.svg)
for (i in 1:2) {
    devs[[i]](width=3, height=2, pointsize=5, file=outs[i])
    layout(1:2)
    par(mar=c(1,4,2,1))
    plot.sbs96(x=0, spectrum=reconstructed.shared.spectrum, main='Shared sSNVs (sigB removed)')
    plot.sbs96(x=0, spectrum=reconstructed.private.spectrum, main='Private sSNVs (sigB removed)')
    dev.off()
}

# raw spectra - i.e., not reconstructed from COSMIC
devs=list(pdf, svglite)
outs=c(out.suppl.pdf, out.suppl.svg)
for (i in 1:2) {
    devs[[i]](width=3, height=2, pointsize=5, file=outs[i])
    layout(1:2)
    par(mar=c(1,4,2,1))
    plot.sbs96(x=0, spectrum=raw.shared.spectrum, main='Shared sSNVs (raw)')
    plot.sbs96(x=0, spectrum=raw.private.spectrum, main='Private sSNVs (raw)')
    dev.off()
}

# barplots of exposure to each signature
devs=list(pdf, svglite)
outs=c(out.expo.pdf, out.expo.svg)
for (i in 1:2) {
    devs[[i]](width=3, height=1.5, pointsize=7, file=outs[i])
    layout(t(1:2))
    par(mar=c(4,4,2,1))
    barplot(sh.expo.pct[names(sh.expo.pct) != 'SigB'],
        main='Shared sSNVs', las=1, border=NA, col='#666666', xlab='Percent of sSNVs', horiz=TRUE)
    barplot(pr.expo.pct[names(sh.expo.pct) != 'SigB'],
        main='Private sSNVs', las=1, border=NA, col='#666666', xlab='Percent of sSNVs', horiz=TRUE)
    dev.off()
}


devs=list(pdf, svglite)
outs=c(out.early_vs_late.pdf)
write.csv(
    cbind(`Early MRCA group`=table(sbs96(muts.early.pairs[shared == TRUE]$mutsig)),
          `Late MRCA group`=table(sbs96(muts.late.pair[shared == TRUE]$mutsig))),
    file=out.early_vs_late.csv, row.names=TRUE, quote=FALSE)
for (i in 1:1) {
    devs[[i]](width=3, height=2.0, pointsize=5, file=outs[i])
    layout(1:2)
    par(mar=c(1,4,2,1))
    plot.sbs96(x=0, spectrum=table(sbs96(muts.early.pairs[shared == TRUE]$mutsig)), main='Early pairs (raw)')
    plot.sbs96(x=0, spectrum=table(sbs96(muts.late.pair[shared == TRUE]$mutsig)), main='Late pair (raw)')
    dev.off()
}

if ('snakemake' %in% ls()) {
    sink()
}
