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
        snakemake@input[1:4], snakemake@output[1:3]
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 7) {
    stop("usage: plot_cosmic_fits.R neuron.expomat.csv neuron.mutburden.csv oligo.expomat.csv oligo.mutburden.csv out.svg out.pdf out.csv")
}

n.expomat.csv <- args[1]
n.mutburden.csv <- args[2]
g.expomat.csv <- args[3]
g.mutburden.csv <- args[4]
out.svg <- args[5]
out.pdf <- args[6]
out.csv <- args[7]

if (file.exists(out.svg))
    stop(paste('output file', out.svg, 'already exists, please delete it first'))
if (file.exists(out.pdf))
    stop(paste('output file', out.pdf, 'already exists, please delete it first'))
if (file.exists(out.csv))
    stop(paste('output file', out.csv, 'already exists, please delete it first'))

suppressMessages(library(data.table))
suppressMessages(library(extrafont))
suppressMessages(library(svglite))

if (!("Arial" %in% fonts()))
    stop("Arial font not detected; did you load extrafonts and run font_import() with the appropriate path?")


nE <- fread(n.expomat.csv)
head(nE)
gE <- fread(g.expomat.csv)
head(gE)
E <- cbind(nE,gE)

nm <- fread(n.mutburden.csv)  # needed for age
gm <- fread(g.mutburden.csv)
m <- rbind(nm, gm)

# Reorder E by increasing age
signames <- E$Sig
E <- E[,-1] # get rid of signature name column

sample.names <- m$sample
E <- E[,..sample.names] # match order of m

n.model.table <-
    do.call(rbind, lapply(1:nrow(E), function(j) {
        data <- data.frame(age=m$age,
            sig.burden=as.matrix(E[j])[1,])[m$type=='PFC neuron',]
        model <- lm(sig.burden ~ age, data=data)
        data.frame(
            CellType='Neuron',
            Sig=signames[j],
            AgeIntercept=unname(coef(model)[1]),
            AgeRate=unname(coef(model)[2]),
            AgeRsquared=summary(model)$r.squared,
            AgeSignif=-log10(coef(summary(model))[2,4]),
            stringsAsFactors=FALSE)
    }))

g.model.table <-
    do.call(rbind, lapply(1:nrow(E), function(j) {
        data <- data.frame(age=m$age,
            sig.burden=as.matrix(E[j])[1,])[m$type=='oligo',]
        model <- lm(sig.burden ~ age, data=data)
        data.frame(
            CellType='Oligo',
            Sig=signames[j],
            AgeIntercept=unname(coef(model)[1]),
            AgeRate=unname(coef(model)[2]),
            AgeRsquared=summary(model)$r.squared,
            AgeSignif=-log10(coef(summary(model))[2,4]),
            stringsAsFactors=FALSE)
    }))


neuron.col <- m$color[m$type=='PFC neuron'][1]
oligo.col <- m$color[m$type=='oligo'][1]

# Handle up to 12 signatures, 3x4 layout
figheight <- 2.5*3
figwidth <- 2.5*4
# Loop over devices to save both pdf and svgs
devs <- list(svglite, pdf)
outs <- c(out.svg, out.pdf)
for (i in 1:2) {
    devs[[i]](file=outs[i], width=figwidth, height=figheight)
    layout(matrix(1:12, nrow=3, byrow=T))
    par(mar=c(4,4,2,1))
    for (j in 1:nrow(E)) {
        plot(m$age, E[j,], pch=17, col=m$color,
            xlab='Age', ylab='Signature exposure', main=signames[j])
        abline(coef=c(n.model.table$AgeIntercept[j], n.model.table$AgeRate[j]),
            lwd=2, col=neuron.col)
        abline(coef=c(g.model.table$AgeIntercept[j], g.model.table$AgeRate[j]),
            lwd=2, col=oligo.col)
        legend('topleft', pch=17, lwd=2, legend=c('Neuron', 'Oligo'),
            col=c(neuron.col, oligo.col))
    }
    dev.off()
}

fwrite(data.table(rbind(n.model.table, g.model.table)), file=out.csv)

if ('snakemake' %in% ls()) {
    sink()
}
