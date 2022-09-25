#!/usr/bin/env Rscript

# detect script being run by snakemake
# if so, make a mock commandArgs function
if ('snakemake' %in% ls()) {
    logfile <- snakemake@log[[1]]
    con <- file(logfile, 'w')
    sink(con, type='output')
    sink(con, type='message')

    # tags come from params
    commandArgs <- function(...) unlist(c(
        snakemake@params[1:2], snakemake@input[1:5], snakemake@output[1]
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 6) {
    cat('celltype and muttype are just text IDs that will be added as 2 columns\n')
    cat('in the output table. e.g., celltype=Neuron, muttype=SNV.\n')
    stop("usage: score_cosmic_signature.R celltype muttype mutmat.csv mutburden.csv cosmic.csv out.csv")
}

tag.celltype <- args[1]
tag.muttype <- args[2]
mutmat.csv <- args[3]
mutburden.csv <- args[4]
cosmic.csv <- args[5]
outcsv <- args[6]

if (file.exists(outcsv))
    stop(paste('output file', outcsv, 'already exists, please delete it first'))

library(pracma)
library(data.table)

cosmic <- fread(cosmic.csv)
cosmic <- cosmic[,-1]  # Get rid of the MutType column (e.g., ACA:C>A)

M <- fread(mutmat.csv)

muttab <- fread(mutburden.csv)

# Remove outliers from M and muttab
muttab <- muttab[outlier == FALSE]
samples.to.keep <- muttab$sample
M <- M[,..samples.to.keep]

###############################################################################
# Evaluate signatures by reduction of fit residual
###############################################################################
# Really simplistic step-forward method for selecting COSMIC signatures.
# Each step selects the sig that minimizes the average residual across
# samples.
# (implying, BTW, that COSMIC fits are being applied per sample.)
#
# IMPORTANT NOTES:
#   1. Each sample is individually fit
#   2. Only minimum residual is computed for each step
top.step=ncol(cosmic)
resids <- c()
terms <- c()
As <- list()  # keep a list of residuals per sample
for (step in 1:top.step) {
    cat(paste('step',step,'\n'))
    print(terms)

    to.test <- (1:ncol(cosmic))
    a <- sapply(to.test, function(j)
        sapply(1:ncol(M), function(i) {
            # important optimization at very late steps (with a huge list of terms)
            if (j %in% terms)
                Inf
            else
                lsqnonneg(as.matrix(cosmic)[,c(j,terms),drop=F], M[[i]])$resid.norm
        }))
    clm <- colMeans(a)
    print(clm)
    s <- which.min(clm)
    print(s)
    terms <- c(terms,s)
    resids <- c(resids, clm[s])
    print(terms)
    print(resids)
}

# The mean residual if no signatures are provided for fitting at all
resid.nosigs <- mean(colSums(M^2))

# sigtable is ordered by the stepwise criteria. that's because
# out of the criteria we compute, only stepwise depends on order.
# start building out the table of metrics for each COSMIC signature
sigtable <- data.table(CellType=tag.celltype,
    MutType=tag.muttype,
    SigName=colnames(cosmic[,..terms]),
    StepwiseMeanResidual=resids,
    StepwiseMeanResidualReduction=-diff(c(resid.nosigs,resids)))

sigtable$PercentStepwiseMeanResidualReduction <-
    100*sigtable$StepwiseMeanResidualReduction / resid.nosigs


###############################################################################
# Evaluate signature stability under sampling (multinomial) noise.
###############################################################################
# Make 10,000 random M matrices (rM)
n.sims=10000
cat('Starting stability simulations:', n.sims, 'random samplings\n')
combined <- rowSums(M)
rM <- rmultinom(n=n.sims, size=sum(combined), prob=combined)
rM.fits <- apply(rM, 2, function(col) lsqnonneg(as.matrix(cosmic), col)$x)
rownames(rM.fits) <- colnames(cosmic)
rM.sd <- apply(rM.fits, 1, sd)
rM.mean <- rowMeans(rM.fits)
rM.cv <- ifelse(rM.mean > 0, rM.sd/rM.mean, 0)  # Coef. of variation

# remember: sigtable is ordered by the stepwise criteria.
# coef. of variation: lower is more stable. use -log10 transform
# to get a higher=more stable metric.
sigtable <- cbind(sigtable,
    StabilityCV=-log10(rM.cv[sigtable$SigName]))


###############################################################################
# Evaluate signatures by total % of mutations contributed.
# E.g., remove signatures contributing <5%, <1%, or whatever cutoff.
###############################################################################
bysample.full <- apply(M, 2, function(col) lsqnonneg(as.matrix(cosmic), col)$x)
rownames(bysample.full) <- colnames(cosmic)
pct.contrib <- 100*rowSums(bysample.full)/sum(bysample.full)
ncells <- rowSums(bysample.full > 0)

sigtable <- cbind(sigtable,
    PercentContribution=pct.contrib[sigtable$SigName],
    NCellsWithSig=ncells[sigtable$SigName])
sigtable$PercentCellsWithSig <- 100*sigtable$NCellsWithSig / ncol(M)


###############################################################################
# Evaluate signatures by correlation with age.
###############################################################################
# Use bysample.full from above. Uses the entire COSMIC catalog.
model.age <- function(burden, age) {
    m <- lm(burden ~ age, data=data.frame(burden=burden, age=age))
    c(AgeRate=unname(coef(m)[2]),
      AgeRsquared=summary(m)$r.squared,
      AgeSignif=-log10(coef(summary(m))[2,4]))
}

mod <- t(apply(bysample.full, 1, model.age, age=muttab$age))

sigtable <- cbind(sigtable, mod[sigtable$SigName,])
sigtable$AgeSignifAdj <- -log10(p.adjust(10^-sigtable$AgeSignif))
                       
fwrite(sigtable, file=outcsv)

if ('snakemake' %in% ls()) {
    sink()
}
