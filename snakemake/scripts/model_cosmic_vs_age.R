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
        as.vector(rbind(snakemake@input[['expomats']],
                        snakemake@input[['mutburdens']]))  # interleave expomats and mutburdens
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 5) {
    stop("usage: plot_cosmic_vs_age.R out.svg out.pdf out.csv expomatrix1.csv mutburden1.csv [expomatrix2.csv mutburden2.csv ... expomatrixN.csv mutburdenN.csv ]")
}

out.svg <- args[1]
out.pdf <- args[2]
out.csv <- args[3]
emat.csvs <- args[seq(4, length(args), 2)]
burden.csvs <- args[seq(5, length(args), 2)]

if (length(emat.csvs) != length(burden.csvs)) {
    stop('each exposure matrix file expomatN.csv must have a matching mutburdenN.csv file')
}

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


signames <- NULL
data <- lapply(1:length(emat.csvs), function(i) {
    E <- fread(emat.csvs[i])

    # get rid of signature name column after saving it once
    if (is.null(signames))
        signames <<- E$Sig
    E <- as.matrix(E[,-1])
    rownames(E) <- signames

    # Reorder E by increasing age
    burden <- fread(burden.csvs[i])[order(age)]  # needed for age
    colnames(burden)[colnames(burden) == 'type'] <- 'celltype'
    sample.names <- burden$sample
    E <- E[,sample.names] # match order of burden
    list(E=E, burden=burden)
})
names(data) <- sapply(data, function(d) d$burden$celltype[1])

combined.data <- list(E=do.call(cbind, lapply(data, function(d) d$E)),
                      burden=do.call(rbind, lapply(data, function(d) d$burden)))

data <- c(data, `All cell types`=list(combined.data))

models <- lapply(data, function(d) {
    rbindlist(lapply(signames, function(signame) {
        data <- data.frame(donor=d$burden$donor, celltype=d$burden$celltype,
                           age=d$burden$age, sig.burden=d$E[signame,])
        celltype <- NA
        if (length(unique(d$burden$celltype)) == 1) {
            model <- lm(sig.burden ~ age, data=data)
            celltype <- d$burden$celltype[1]
        } else {
            model <- lm(sig.burden ~ age*celltype, data=data)
            celltype <- 'All cell types'
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
            CellType=celltype,
            Sig=signame,
            Formula=deparse(formula(model)),
            Variable=c('AgeRsquared', rownames(coefs)),
            # Add Rsquared as a variable
            rbind(c(summary(model)$r.squared, rep(NA, ncol(coefs)+2-1)),
                  cbind(coefs, ci[rownames(coefs),]))
        )
    }))
})

fwrite(rbindlist(models), file=out.csv) #data.table(rbind(n.model.table, g.model.table)), file=out.csv)

# Handle up to 16 signatures, 4x4 layout
figheight <- 2.5*4
figwidth <- 2.5*4
# Loop over devices to save both pdf and svgs
devs <- list(svglite, pdf)
outs <- c(out.svg, out.pdf)
for (i in 1:2) {
    devs[[i]](file=outs[i], width=figwidth, height=figheight)
    layout(matrix(1:16, nrow=4, byrow=T))
    par(mar=c(4,4,2,1))
    for (signame in signames) {
        print(signame)
        all.data <- data[[length(data)]]
        plot(all.data$burden$age, all.data$E[signame,], pch=17,
            col=all.data$burden$color,
            xlab='Age', ylab='Signature exposure', main=signame)
        for (i in 1:(length(models) - 1)) { # don't plot the combined model
            m <- models[[i]][Sig == signame]
            abline(coef=c(m[Variable == '(Intercept)']$Estimate,
                          m[Variable == 'age']$Estimate),
                lwd=2, col=data[[i]]$burden$color[1])
        }
        legend('topleft', pch=17, lwd=2, legend=names(models)[-length(models)],
            col=sapply(data[-length(data)], function(d) d$burden$color[1]))
    }
    dev.off()
}


if ('snakemake' %in% ls()) {
    sink()
}
