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
        paste(snakemake@params[['colors']], collapse=','),
        paste(unlist(snakemake@params['group_names']), snakemake@input, sep='=')
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 6) {
    cat('`colors` is a comma separated string of R recognized color strings IN THE SAME ORDER AS THE burden_files.\n')
    cat("`colors=table` to use the colors already present in the burden files.\n")
    stop("usage: model_aging_rates.r out.svg out.pdf out.combined_burdens.csv out.models.csv colors group1_name=burden_file1.csv [ group2_name=burden_file2.csv ... groupN_name=burden_fileN.csv ]")
}

outsvg <- args[1]
outpdf <- args[2]
outcsv <- args[3]
outmodelcsv <- args[4]
color.string <- args[5]
burden.csvs <- args[-(1:5)]

if (color.string != "table") {
    colors <- unlist(strsplit(color.string, split=',')[[1]])
    print(colors)
} else {
    colors <- color.string
    print("colors=table, not overriding table colors")
}

if (file.exists(outsvg))
    stop(paste('output file', outsvg, 'already exists, please delete it first'))
if (file.exists(outpdf))
    stop(paste('output file', outpdf, 'already exists, please delete it first'))
if (file.exists(outcsv))
    stop(paste('output file', outcsv, 'already exists, please delete it first'))
if (file.exists(outmodelcsv))
    stop(paste('output file', outmodelcsv, 'already exists, please delete it first'))

group.names <- sapply(strsplit(burden.csvs, '='), head, 1)
burden.csvs <- sapply(strsplit(burden.csvs, '='), function(x) paste(x[-1], collapse='='))

suppressMessages(library(data.table))
suppressMessages(library(lme4))
suppressMessages(library(lmerTest))
suppressMessages(library(svglite))


burdens <- setNames(lapply(1:length(burden.csvs), function(i) {
    in.csv <- burden.csvs[i]
    burden <- fread(in.csv)
    burden[, group := group.names[i]]
    if (colors != "table") {
        cat(paste0('group=', group.names[i], ': overiding table color=', burden$color[1], ' with color=', colors[i], '\n'))
        burden[, color := colors[i]]
    }
    if (!('uncorrected.genome.burden' %in% colnames(burden))) {
        cat(paste0('group=', group.names[i], ': adding copyover uncorrected.genome.burden column to table'))
        burden[, uncorrected.genome.burden := genome.burden]
    }
    burden
}), group.names)
print(burdens)

# harmonize column names/order before rbindlist
# unique maintains order
all.col.names <- unique(unlist(lapply(burdens, colnames)))
burdens <- lapply(burdens, function(b) setcolorder(b, all.col.names))

combined.burdens <- rbindlist(burdens)
print(combined.burdens)

burdens <- c(burdens, `All groups`=list(combined.burdens))
print(burdens)
models <- lapply(burdens, function(dt) {
    if (length(unique(dt$group)) == 1) {
        model <- lmer(genome.burden ~ age + (1|donor), data=dt[outlier == "NORMAL"])
    } else {
        # for the combined model
        model <- lmer(genome.burden ~ age*group + (1|donor), data=dt[outlier == "NORMAL"])
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
        predict(m$model, data.frame(age=max(combined.burdens$age), donor=m$model@frame$donor[1]))
    )
    # make ylim big enough to contain all data points and lines
    ylim <- c(0, max(combined.burdens$genome.burden, maxpred, na.rm=TRUE))
    plot(combined.burdens$plotage,
        combined.burdens$genome.burden,
        col=combined.burdens$color, pch=ifelse(combined.burdens$outlier == 'NORMAL', 17, 4),
        ylim=ylim, ylab='Autosomal mutation burden', xlab='Age')
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
    data.table(Group=names(models)[i],
        Formula=deparse(formula(m$model)),
        Variable=rownames(coefs),
        cbind(coefs, m$ci[rownames(coefs),]))
}))
fwrite(all.models, file=outmodelcsv)

if ('snakemake' %in% ls()) {
    sink()
}
