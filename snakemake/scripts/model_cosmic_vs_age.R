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
        snakemake@output['svg'],
        snakemake@output['pdf'],
        snakemake@output['csv'],
        snakemake@input['metadata'],
        paste(snakemake@params[['colors']], collapse=','),
        paste(snakemake@params[['groups']], snakemake@input[['expomats']], sep='=')
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 5) {
    cat("colors is a comma separated string of R recognized colors. Order must match expomatrix CSV files.\n")
    cat("`colors=table` to use the colors already present in the metadata file.\n")
    stop("usage: model_cosmic_vs_age.R out.svg out.pdf out.csv metadata.csv colors expomatrix1.csv [ expomatrix2.csv ... expomatrixN.csv ]")
}

out.svg <- args[1]
out.pdf <- args[2]
out.csv <- args[3]
metadata.csv <- args[4]
color.string <- args[5]
emat.csvs <- args[-(1:5)]

if (file.exists(out.svg))
    stop(paste('output file', out.svg, 'already exists, please delete it first'))
if (file.exists(out.pdf))
    stop(paste('output file', out.pdf, 'already exists, please delete it first'))
if (file.exists(out.csv))
    stop(paste('output file', out.csv, 'already exists, please delete it first'))

if (color.string != "table") {
    colors <- unlist(strsplit(color.string, split=',')[[1]])
    print(colors)
} else {
    colors <- color.string
    print("colors=table, not overriding table colors")
}

suppressMessages(library(data.table))
suppressMessages(library(svglite))

# Group names are only relevant for the combined model
group.names <- sapply(strsplit(emat.csvs, '='), head, 1)
emat.csvs <- sapply(strsplit(emat.csvs, '='), function(x) paste(x[-1], collapse='='))

meta <- fread(metadata.csv)
meta[, group := '']

signames <- NULL
data <- setNames(lapply(1:length(emat.csvs), function(i) {
    E <- fread(emat.csvs[i])

    # get rid of signature name column after saving it once
    # implies all groups must be fit to the same set of signatures.
    if (is.null(signames))
        signames <<- E$Sig
    E <- as.matrix(E[,-1])
    rownames(E) <- signames

    sample.names <- meta[order(age)][sample %in% colnames(E)]$sample
    meta[sample %in% sample.names, group := group.names[i]]
    if (colors != 'table')
        meta[sample %in% sample.names, color := colors[i]]
    E
}), group.names)
combined.data <- do.call(cbind, data)
data <- c(data, `All groups`=list(combined.data))

# Assumes each group has only one color in it
if (colors == 'table') {
    group.colors <- setNames(sapply(group.names, function(name) meta[group == name]$color[1]), group.names)
} else {
    group.colors <- setNames(colors, group.names)
}
print(group.colors)


models <- lapply(data, function(E) {
    rbindlist(lapply(signames, function(signame) {
        data <- meta[data.table(sample=colnames(E), sig.burden=E[signame,]),,on=.(sample)]
        group <- NA
        if (length(unique(data$group)) == 1) {
            model <- lm(sig.burden ~ age, data=data)
            group <- data$group[1]
        } else {
            model <- lm(sig.burden ~ age*group, data=data)
            group <- 'All groups'
        }

        ci <- confint(model)
        colnames(ci) <- paste0('95% CI ', c('lower', 'upper'))

        coefs <- coef(summary(model))
        # this column name is annoying to access programmatically - change it
        colnames(coefs)[colnames(coefs) == 'Pr(>|t|)'] <- 'Pvalue'
        # only adjust P-values for slopes - we don't make any claims about intercepts
        coefs <- cbind(coefs, Padj=NA)
        slope.rows <- which(substr(rownames(coefs), 1, 3) == 'age')
        coefs[slope.rows, 'Padj'] <- p.adjust(coefs[slope.rows, 'Pvalue'])
        data.table(
            Group=group,
            Sig=signame,
            Formula=deparse(formula(model)),
            Variable=c('AgeRsquared', rownames(coefs)),
            # Add Rsquared as a variable
            rbind(c(summary(model)$r.squared, rep(NA, ncol(coefs)+2-1)),
                  cbind(coefs, ci[rownames(coefs),]))
        )
    }))
})

fwrite(rbindlist(models), file=out.csv)

# Handle up to 20 signatures, 4x5 layout
figheight <- 2.5*4
figwidth <- 2.5*5
# Loop over devices to save both pdf and svgs
devs <- list(svglite, pdf)
outs <- c(out.svg, out.pdf)
for (i in 1:2) {
    devs[[i]](file=outs[i], width=figwidth, height=figheight)
    layout(matrix(1:20, nrow=4, byrow=T))
    par(mar=c(4,4,2,1))
    for (signame in signames) {
        print(signame)
        all.E <- data[[length(data)]] # Get the exposure matrix with all samples present
        all.data <- meta[data.table(sample=colnames(all.E), sig.burden=all.E[signame,]),,on=.(sample)]
        plot(all.data$age, all.data$sig.burden,
            pch=ifelse(all.data$outlier == 'NORMAL', 17, 4),
            col=all.data$color,
            xlab='Age', ylab='Signature exposure', main=signame)
        for (i in 1:(length(models) - 1)) { # don't plot the combined model
            m <- models[[i]][Sig == signame]
            abline(coef=c(m[Variable == '(Intercept)']$Estimate,
                          m[Variable == 'age']$Estimate),
                lwd=2, col=group.colors[group.names[i]])
        }
        legend('topleft', pch=17, lwd=2, legend=names(models)[-length(models)],
            col=group.colors[names(models)])
    }
    dev.off()
}


if ('snakemake' %in% ls()) {
    sink()
}
