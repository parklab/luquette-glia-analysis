#!/usr/bin/env Rscript

# detect script being run by snakemake
# if so, make a mock commandArgs function
if ('snakemake' %in% ls()) {
    logfile <- snakemake@log[[1]]
    con <- file(logfile, 'w')
    sink(con, type='output')
    sink(con, type='message')

    commandArgs <- function(...) unlist(c(
        snakemake@params['muttype'],
        snakemake@params['amptype'],
        snakemake@input['muts'],
        snakemake@input['mutburden'],
        snakemake@input['cosmic'],
        snakemake@input['sigb'],
        snakemake@output['mutmat'],
        snakemake@output['expomat']
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 8) {
    cat('muttype must be either SNV or Indel\n')
    cat('when amptype=MDA (case sensitive), Lodato et al. signature B is added to COSMIC, but its exposure is REMOVED from the output\n')
    stop("usage: analyze_cosmic_aging.R muttype amptype mutations.csv mutburden.csv cosmic.csv out.mutmat.csv out.expo.csv")
}

muttype <- args[1]
amptype <- args[2]
inmuts <- args[3]
inmutburden <- args[4]
incosmic <- args[5]
insigb <- args[6]
out.mutmat.csv <- args[7]
out.expo.csv <- args[8]

if (muttype != 'SNV' & muttype != 'Indel')
    stop(paste('muttype must be either SNV or Indel, got', muttype))
if (amptype != 'PTA' & amptype != 'MDA')
    stop(paste('amptype must be either PTA or MDA (case sensitive), got', amptype))
if (file.exists(out.mutmat.csv))
    stop(paste('output file', out.mutmat.csv, 'already exists, please delete it first'))
if (file.exists(out.expo.csv))
    stop(paste('output file', out.expo.csv, 'already exists, please delete it first'))

library(scan2)       # for df.to.sbs96
library(pracma)

# Given a list of mutation dataframes, calculate exposure to the supplied
# COSMIC database (which may have been subsetted to contain only spectra
# passing stability thresholds).
# COSMIC format: column 1 is the mutation type. Need to remove it before
# fitting.
exposure <- function(x, cosmic=cosmic) {
    # lsqnonneg fails if given NA
    if (any(is.na(x)))
        ret <- rep(NA, ncol(cosmic))
    else
        ret <- lsqnonneg(as.matrix(cosmic), x)$x
    setNames(ret, colnames(cosmic))
}


muts <- fread(inmuts)
mutburden <- fread(inmutburden)  # allows extrapolation to genome-wide numbers
cosmic <- fread(incosmic)
sigb <- fread(insigb)$B

if (amptype == 'MDA') {
    cat("amptype=MDA, adding Lodato et al. Signature B to COSMIC catalog\n")
    cosmic$SigB <- sigb
}

# nsom=0 will unfortunately get extrapolated to 0
# genome.burden=NA for outlier=TRUE samples. need to propagate NA
# so that the measurement is not considered by linear models but
# cannot allow certain functions (lsqnonneg) to see NA (they fail)
mutburden[, correction.factor := ifelse(nsom == 0, 0, genome.burden / nsom)]

# Get the finalized mutation matrix M, which involves mapping mutations to
# SBS96 or ID83 channels and multiplying by the genome-wide burden correction factor.
# M is [96|83] x (1 + #samples), sample order is the same as mutburden
samples <- mutburden[['sample']]
if (muttype == 'SNV') {
    M <- sapply(samples, function(s)
        table(sbs96(muts[sample==s,]$mutsig)))
} else {
    # plot.indel returns the ID83 signature vector
    M <- sapply(samples, function(s)
        table(id83(muts[sample==s,]$mutsig)))
}
head(M)

M <- M[,mutburden[['sample']]]  # reorder to match mutburden
M <- t(t(M) * mutburden[['correction.factor']])

# sanity check to ensure matching orders
if (!all(rownames(M) == cosmic[,1]))
    stop('mutations and COSMIC signatures are not in the same order')

str(M)

head(cosmic)
str(as.matrix(cosmic[,-1]))
# make sure not to pass mutation type columns
E <- apply(M, 2, exposure, cosmic=as.matrix(cosmic[,-1]))
# data.tables don't use rownames
E <- data.table(Sig=rownames(E), E)
head(E)

if (amptype == 'MDA') {
    cat("Removing Signature B contribution from exposure matrix\n")
    E <- E[Sig != 'SigB',]
}

fwrite(M, file=out.mutmat.csv)
fwrite(E, file=out.expo.csv)

if ('snakemake' %in% ls()) {
    sink()
}
