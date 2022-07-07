#!/usr/bin/env Rscript

# detect script being run by snakemake
# if so, make a mock commandArgs function
if ('snakemake' %in% ls()) {
    logfile <- snakemake@log[[1]]
    con <- file(logfile, 'w')
    sink(con, type='output')
    sink(con, type='message')

    commandArgs <- function(...) unlist(c(
        snakemake@input[1:2], snakemake@output[1:3], snakemake@params
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 12) {
    stop("usage: rda_to_csv.r neuron_scores.csv oligo_scores.csv out.svg out.pdf out.csv pct_resid_reduc_cutoff stability_cutoff pct_contrib_cutoff pct_cells_cutoff muts_per_year_cutoff age_signif_cutoff num_cutoffs_cutoff [ sig(s)_to_exclude ... ]")
}

nscore.csv <- args[1]
gscore.csv <- args[2]
out.svg <- args[3]
out.pdf <- args[4]
out.csv <- args[5]
resid.reduc.cutoff <- as.numeric(args[6])
stab.cutoff <- as.numeric(args[7])
pct.contrib.cutoff <- as.numeric(args[8])
pct.cells.cutoff <- as.numeric(args[9])
muts.per.year.cutoff <- as.numeric(args[10])
age.signif.cutoff <- as.numeric(args[11])
numcutoffs.cutoff <- as.numeric(args[12])
sigs.to.exclude <- args[-(1:12)]

if (file.exists(out.svg))
    stop(paste('outut file', out.svg, 'already exists, please delete it first'))
if (file.exists(out.pdf))
    stop(paste('outut file', out.pdf, 'already exists, please delete it first'))
if (file.exists(out.csv))
    stop(paste('outut file', out.csv, 'already exists, please delete it first'))

suppressMessages(library(data.table))
suppressMessages(library(extrafont))
suppressMessages(library(svglite))

if (!("Arial" %in% fonts()))
    stop("Arial font not detected; did you load extrafonts and run font_import() with the
 appropriate path?")

n <- fread(nscore.csv)
g <- fread(gscore.csv)

# Make combined table for taking maximums. The neuron and oligo tables may
# not be in the same order, since they're eached ordred by the stepwise
# residual reduction method.
ntmp <- copy(n[,-(1:2)])  # columns 1+2 are celltype and muttype
setkey(ntmp, SigName)
colnames(ntmp)[-1] <- paste0('Neuron', colnames(ntmp)[-1])
gtmp <- copy(g[,-(1:2)])
setkey(gtmp, SigName)
colnames(gtmp)[-1] <- paste0('Oligo', colnames(gtmp)[-1])
combined <- ntmp[gtmp]  # Join
combined[, `:=`(
        PercentStepwiseMeanResidualReduction = pmax(NeuronPercentStepwiseMeanResidualReduction, OligoPercentStepwiseMeanResidualReduction),
        PercentCellsWithSig = pmax(NeuronPercentCellsWithSig, OligoPercentCellsWithSig),
        StabilityCV = pmax(NeuronStabilityCV, OligoStabilityCV),
        PercentContribution = pmax(NeuronPercentContribution, OligoPercentContribution),
        AgeRate = pmax(NeuronAgeRate, OligoAgeRate),
        AgeSignif = pmax(NeuronAgeSignif, OligoAgeSignif))]
# Add cutoffs
combined[, `:=`(
        CutoffPercentStepwiseMeanResidualReduction = PercentStepwiseMeanResidualReduction >= resid.reduc.cutoff,
        CutoffStabilityCV = StabilityCV >= stab.cutoff,
        CutoffPercentContribution = PercentContribution >= pct.contrib.cutoff,
        CutoffPercentCellsWithSig = PercentCellsWithSig >= pct.cells.cutoff,
        CutoffAgeRate = AgeRate >= muts.per.year.cutoff,
        CutoffAgeSignif = AgeSignif >= age.signif.cutoff,
        SigNamedExclusion = SigName %in% sigs.to.exclude)]
combined[, NumCutoffs := CutoffPercentStepwiseMeanResidualReduction +
                         CutoffStabilityCV +
                         CutoffPercentContribution +
                         CutoffPercentCellsWithSig +
                         CutoffAgeRate +
                         CutoffAgeSignif]
combined[, SigIncluded := NumCutoffs >= numcutoffs.cutoff & !SigNamedExclusion]
head(combined)
table(combined$NumCutoffs)
combined[SigIncluded == TRUE]

# Make the set of histograms to aid cutoff selection
plot6 <- function(d, border=NA, col='#4c4c4c', breaks=20, ylab='Number of signatures', ...) {
    par(mar=c(4,4,2,1))
    hist(log10(d$PercentStepwiseMeanResidualReduction),
        breaks=breaks, ylab=ylab, col=col, border=NA,
        xlab='log10(Percent residual reduction)', ...)
    abline(v=log10(resid.reduc.cutoff), col=2, lwd=2)
    hist(d$StabilityCV, breaks=breaks, ylab=ylab, col=col, border=NA,
        xlab='Stability (multinomial noise)', ...)
    abline(v=stab.cutoff, col=2, lwd=2)
    hist(d$PercentContribution, breaks=breaks, ylab=ylab, col=col, border=NA,
        xlab='Percent contribution', ...)
    abline(v=pct.contrib.cutoff, col=2, lwd=2)
    hist(d$PercentCellsWithSig, breaks=breaks, ylab=ylab, col=col, border=NA,
        xlab='Percent cells with sig.', ...)
    abline(v=pct.cells.cutoff, col=2, lwd=2)
    hist(d$AgeRate, breaks=breaks, ylab=ylab, col=col, border=NA,
        xlab='Mutations per year', ...)
    abline(v=muts.per.year.cutoff, col=2, lwd=2)
    hist(d$AgeSignif, breaks=breaks, ylab=ylab, col=col, border=NA,
        xlab='-log10(P-value in age model)', ...)
    abline(v=age.signif.cutoff, col=2, lwd=2)
}


figwidth=14
figheight=9
devs <- list(svglite, pdf)
outs <- c(out.svg, out.pdf)
for (i in 1:2) {
    cat('writing', outs[i], '\n')
    devs[[i]](file=outs[i], width=figwidth, height=figheight)
    layout(matrix(1:(3*6),nrow=3,byrow=TRUE))
    plot6(n, main='Neurons')
    plot6(g, main='Oligo')
    plot6(combined, main='Max(Neurons, Oligo)')

    dev.off()
}

fwrite(combined, file=out.csv)

if ('snakemake' %in% ls()) {
    sink()
}
