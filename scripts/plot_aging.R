#!/usr/bin/env Rscript

# detect script being run by snakemake
# if so, make a mock commandArgs function
if ('snakemake' %in% ls()) {
    logfile <- snakemake@log[[1]]
    con <- file(logfile, 'w')
    sink(con, type='output')
    sink(con, type='message')

    commandArgs <- function(...) unlist(c(
        snakemake@input[1:2], snakemake@output[1:4]
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 6) {
    cat('assumes in.rda contains only one object\n')
    stop("usage: rda_to_csv.r neuron_burden.csv oligo_burden.csv out.svg out.pdf out_table.csv out_model.csv")
}

nfile <- args[1]
gfile <- args[2]
outsvg <- args[3]
outpdf <- args[4]
outcsv <- args[5]
outmodelcsv <- args[6]

if (file.exists(outsvg))
    stop(paste('output file', outsvg, 'already exists, please delete it first'))
if (file.exists(outpdf))
    stop(paste('output file', outpdf, 'already exists, please delete it first'))
if (file.exists(outcsv))
    stop(paste('output file', outcsv, 'already exists, please delete it first'))
if (file.exists(outmodelcsv))
    stop(paste('output file', outmodelcsv, 'already exists, please delete it first'))

suppressMessages(library(data.table))
suppressMessages(library(lme4))
suppressMessages(library(lmerTest))
suppressMessages(library(extrafont))
suppressMessages(library(svglite))

if (!("Arial" %in% fonts()))
    stop("Arial font not detected; did you load extrafonts and run font_import() with the appropriate path?")


n <- fread(nfile)
g <- fread(gfile)
combined <- rbind(n, g)

nmodel <- lmer(genome.burden ~ age + (1|donor), data=n)
gmodel <- lmer(genome.burden ~ age + (1|donor), data=g)


figwidth=4
figheight=4
devs <- list(svglite, pdf)
outs <- c(outsvg, outpdf)
for (i in 1:2) {
    devs[[i]](file=outs[i], width=figwidth, height=figheight)

    # any donor can be used
    nline <- predict(nmodel, data.frame(age=max(combined$age), donor=combined$donor[1]))
    gline <- predict(gmodel, data.frame(age=max(combined$age), donor=combined$donor[1]))
    ylim <- c(0, max(nline, gline, combined$genome.burden))
    plot(combined$plotage,
        combined$genome.burden,
        col=combined$color, pch=17, ylim=ylim,
        ylab='Autosomal mutation burden', xlab='Age')
    abline(coef=coef(summary(nmodel)), lwd=2, col=n$color[1])
    abline(coef=coef(summary(gmodel)), lwd=2, col=g$color[2])
    legend('topleft', lwd=2, pch=17, col=1:2, legend=c('Neurons', 'Oligo'))
}

fwrite(combined, file=outcsv)
fwrite(rbind(
    data.table(CellType='Neuron', Variable=c('Intercept','Age'), coef(summary(nmodel))),
    data.table(CellType='Oligo', Variable=c('Intercept','Age'), coef(summary(gmodel)))
), file=outmodelcsv)


if ('snakemake' %in% ls()) {
    sink()
}
