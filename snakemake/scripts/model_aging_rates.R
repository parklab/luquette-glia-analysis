#!/usr/bin/env Rscript

# detect script being run by snakemake
# if so, make a mock commandArgs function
if ('snakemake' %in% ls()) {
    logfile <- snakemake@log[[1]]
    con <- file(logfile, 'w')
    sink(con, type='output')
    sink(con, type='message')

    commandArgs <- function(...) unlist(c(
        snakemake@output['svg'],
        snakemake@output['pdf'],
        snakemake@output['burdens'],
        snakemake@output['models'],
        snakemake@input
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 5) {
    cat("IMPORTANT: burden_fileX.csv must contain a `type` column that specifies the cell type of each sample.")
    stop("usage: model_aging_rates.r out.svg out.pdf out.combined_burdens.csv out.models.csv burden_file1.csv [ burden_file2.csv ... burden_fileN.csv ]")
}

outsvg <- args[1]
outpdf <- args[2]
outcsv <- args[3]
outmodelcsv <- args[4]
burden.csvs <- args[-(1:4)]

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
suppressMessages(library(svglite))


burdens <- lapply(burden.csvs, function(in.csv) {
    burden <- fread(in.csv)
    if (!('type' %in% colnames(burden))) {
        stop(paste0('burden_file (', in.csv, ') must contain a `type` column'))
    }
    colnames(burden)[colnames(burden) == 'type'] <- 'celltype'
    if (length(unique(burden$celltype)) > 1) {
        stop(paste0('all lines in burden_file (', in.csv, ') must contain the same `type` value'))
    }
    burden
})
names(burdens) <- sapply(burdens, function(b) b$celltype[1])
print(burdens)

combined.burdens <- rbindlist(burdens)
print(combined.burdens)

burdens <- c(burdens, `All cell types`=list(combined.burdens))
print(burdens)
models <- lapply(burdens, function(dt) {
    if (length(unique(dt$celltype)) == 1) {
        model <- lmer(genome.burden ~ age + (1|donor), data=dt)
    } else {
        # for the combined model
        model <- lmer(genome.burden ~ age*celltype + (1|donor), data=dt)
    }
    ci <- confint(model)
    colnames(ci) <- paste0('95% CI ', c('lower', 'upper'))
    list(model=model, ci=ci)
})
print(models)

figwidth=4
figheight=4
devs <- list(svglite, pdf)
outs <- c(outsvg, outpdf)
for (i in 1:2) {
    devs[[i]](file=outs[i], width=figwidth, height=figheight)

    # any donor can be used
    # [-length(models)] - don't use the combined model
    maxpred <- sapply(models[-length(models)], function(m)
        predict(m$model, data.frame(age=max(combined.burdens[outlier==FALSE]$age), donor=m$model@frame$donor[1]))
    )
    # make ylim big enough to contain all data points and lines
    ylim <- c(0, max(combined.burdens$genome.burden, maxpred, na.rm=TRUE))
    plot(combined.burdens$plotage,
        combined.burdens$genome.burden,
        col=combined.burdens$color, pch=17, ylim=ylim,
        ylab='Autosomal mutation burden', xlab='Age')
    abline(h=0, col='grey')
    # These curve() calls don't plot below 0.
    for (i in 1:(length(models)-1)) {
        m <- models[[i]]
        curve(coef(summary(m$model))[1] + coef(summary(m$model))[2]*x,
            from=0, to=2*max(combined.burdens$age), lwd=2,
            col=burdens[[i]]$color[1], add=T)
    }
    legend('topleft', lwd=2, pch=17,
        col=sapply(burdens[-length(burdens)], function(b) b$color[1]),
        legend=names(burdens)[-length(burdens)])
}

fwrite(combined.burdens, file=outcsv)

for (i in 1:length(models)) {
    m <- models[[i]]
    print(names(models)[i])
    print(summary(m$model))
}

all.models <- rbindlist(lapply(1:length(models), function(i) {
    m <- models[[i]]
    coefs <- coef(summary(m$model))
    data.table(CellType=names(models)[i],
        Formula=deparse(formula(m$model)),
        Variable=rownames(coefs),
        cbind(coefs, m$ci[rownames(coefs),]))
}))
fwrite(all.models, file=outmodelcsv)

if ('snakemake' %in% ls()) {
    sink()
}
