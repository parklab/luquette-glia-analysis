#!/usr/bin/env Rscript

# detect script being run by snakemake
# if so, make a mock commandArgs function
if ('snakemake' %in% ls()) {
    logfile <- snakemake@log[[1]]
    con <- file(logfile, 'w')
    sink(con, type='output')
    sink(con, type='message')

    commandArgs <- function(...) {
        ret <- unlist(c(
            snakemake@input[1],         # table from make_roadmap_enrichment_table.R
            snakemake@params['ignore'],
            snakemake@output[1:3]      # output svg/pdf/csv ||no longer: /jpeg
        ))
        ret
    }
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 2) {
    cat("enrichment_table.csv is expected to have metadata columns in the order:\n")
    cat("  1. QUANTILES, 2. BINSIZE, 3-(n), metadata, (n+1)-end, esummary() rows.\n")
    cat("Any number of metadata columns is allowed; metadata column names cannot\n")
    cat("collide with esummary() column names.\n")
    stop("usage: make_roadmap_enrichment_grid.R out.csv enrichment_table1.csv [ enrichment_table2.csv ... enrichment_tableN.csv ]")
}

out.csv <- args[1]
in.csvs<- args[-1]

if (file.exists(out.csv))
    stop(paste('output file', out.csv, 'already exists, please delete it first'))

library(data.table)
library(mutenrich)

ds <- lapply(in.csvs, function(in.csv) {
    cat(in.csv)
    d <- fread(in.csv)

    d$QUANTILES <- as.integer(d$QUANTILES)
    d$BINSIZE <- as.integer(d$BINSIZE)

    # Get rid of the 'excluded' and 'outside' classes, which have non-integer quantile values
    d <- d[!is.na(as.integer(d$quantile))]
    d$quantile <- as.integer(d$quantile)
    
    # use nr. of quantiles and binsize as metadata now
    meta.cols <- colnames(d)
    meta.cols <- meta.cols[(1:length(meta.cols)) < which(meta.cols=='quantile')]
    # if "sigclass" is present (from sigenrich analysis), use it as the
    # primary data grouping factor.
    if ('sigclass' %in% meta.cols)
        meta.cols <- c('sigclass', meta.cols[-which(meta.cols=='sigclass')])
    group.factors <- paste(meta.cols, collapse=',')
    groups <- list(no.name=d)
    if (length(group.factors) > 0) {
        group.factors <- strsplit(group.factors, ',')[[1]]
        # Determine the columns by group.factors
        cat('group.factors', group.factors, '\n')
        for (gf in group.factors) {
            if (!(gf %in% colnames(d)))
                stop(paste('group.factor', gf, 'is not a valid column name.\nAvailable column names:', colnames(d)))
        }
        groups <- split(d, d[, apply(.SD, 1, paste, collapse=','), .SDcols=group.factors])
    }

    cat("Got", length(groups), "groups (each = 1 plot column):\n")

    do.call(rbind, lapply(1:length(groups), function(idx) {
        group.name <- names(groups)[idx]

        g <- groups[[group.name]]
        # recompute this to remove the quantile and binsize info from the meta tag
        meta.cols <- colnames(g)[!(colnames(g) %in% c('QUANTILES', 'BINSIZE'))]
        meta.cols <- meta.cols[(1:length(meta.cols)) < which(meta.cols=='quantile')]
        meta <- paste(g[1,..meta.cols], collapse=',')
        npoints <- sum(!is.na(g$enr))
        scores <- c(NA, NA)
        if (npoints > 1) {  # need at least 2 points for regression. some covariates are too sparse and have no different quantiles
            m <- lm(enr ~ quantile, data=g)
#print(g)
#print(summary(m))
            # return both the slope and the p-value of the vs. 0 test
            scores <- coef(summary(m))[2,c(1,4)]
        }

        ret <- data.frame(nquantiles=g$QUANTILES[1], binsize=g$BINSIZE[1], npoints=npoints,
            meta=meta, slope=scores[1], pval=scores[2])
        return(ret)
    }))

})

fwrite(data.table(do.call(rbind, ds)), file=out.csv)

if ('snakemake' %in% ls()) {
    sink()
}
