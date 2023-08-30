#!/usr/bin/env Rscript

# detect script being run by snakemake
# if so, make a mock commandArgs function
if ('snakemake' %in% ls()) {
    logfile <- snakemake@log[[1]]
    con <- file(logfile, 'w')
    sink(con, type='output')
    sink(con, type='message')

    commandArgs <- function(...) unlist(c(
        snakemake@input['cosmic'],
        snakemake@output['csv'],
        snakemake@params['add_pta']
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 3) {
    cat("add_pta: add the PTA artifact signature\n")
    stop("usage: normalize_cosmic_snv_table_and_add_pta_artifact.R in_cosmic.csv out.csv {add_pta: TRUE|FALSE}")
}


cosmic.csv <- args[1]
outcsv <- args[2]
add.pta <- as.logical(args[3])

if (file.exists(outcsv))
    stop(paste('output file', outcsv, 'already exists, please delete it first'))

suppressMessages(library(scan2))

cosmic <- fread(cosmic.csv, header=T, stringsAsFactors=F)
head(cosmic)
if (colnames(cosmic)[1] == 'MutationType') {
    # Format of signature active signature catalog estimated by SigProfilerExtractor
    # Convert SigProfilerExtractor format into the format of signatures downloaded
    # from the COSMIC website.
    cat("Converting SigProfilerExtractor output to COSMIC format..\n")
    cosmic <- cbind(Type=substr(cosmic[[1]], 3, 5),
                    Subtype=paste0(substr(cosmic[[1]], 1, 1), substr(cosmic[[1]], 3, 3), substr(cosmic[[1]], 7, 7)),
                    cosmic[,-1])
    print(head(cosmic))
}

# older code used read.csv instead of fread. just reverting to data.frame to allow
# things to run unaltered.
cosmic <- as.data.frame(cosmic)

rownames(cosmic) <- paste(cosmic[,2], cosmic[,1], sep=':')

# These are the columns that have " -    " values in them.
str(cosmic)
#print(cosmic$SBS25)
#print(cosmic$SBS27)
#print(cosmic$SBS46)
#print(cosmic$SBS47)
#print(cosmic$SBS84)
#print(cosmic$SBS88)
# A few columns have " -   " values for some trinuc contexts. Not sure why.
# By looking at these signatures on the COSMIC website, it seems these
# values just mean 0.
#
# UPDATE: now using fread(), which trims the extra spaces around -
cosmic[cosmic == "-"] <- 0
tmp.rownames <- rownames(cosmic)
cosmic <- as.data.frame(lapply(cosmic[-(1:2)], function(column) as.numeric(column)),
    stringsAsFactors=FALSE)
rownames(cosmic) <- tmp.rownames
head(cosmic)

# reorder SBS channels to match the plot order (which is also SCAN2's order)
o <- as.character(names(table(sbs96(c()))))
cosmic <- cosmic[o,]
cosmic <- cbind(MutType=rownames(cosmic), cosmic)

if (add.pta) {
    # Add the PTA error signature
    lysis.sig <- t(t(get(data(snv.artifact.signature.v3))))

    if (!all(rownames(lysis.sig) == rownames(cosmic))) {
        cat('pta order:'); print(rownames(lysis.sig))
        cat('cosmic order:'); print(rownames(cosmic))
        stop("PTA artifact signature and COSMIC mutation types are not in the same order")
    }
    cosmic <- cbind(cosmic, lysis.sig)
}


fwrite(cosmic, file=outcsv)

if ('snakemake' %in% ls()) {
    sink()
}
