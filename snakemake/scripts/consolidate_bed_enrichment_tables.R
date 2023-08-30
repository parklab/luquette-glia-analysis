#!/usr/bin/env Rscript

# detect script being run by snakemake
# if so, make a mock commandArgs function
if ('snakemake' %in% ls()) {
    logfile <- snakemake@log[[1]]
    con <- file(logfile, 'w')
    sink(con, type='output')
    sink(con, type='message')

    commandArgs <- function(...) unlist(c(
        snakemake@params['group_tag'],
        snakemake@params['signal'],
        snakemake@params['features'],
        snakemake@output['csv'],
        snakemake@input['enrich_objects']
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 5) {
    cat("tag: arbitrary string to identify the group of cells these enrichment objects were derived from\n")
    cat('signal: specified as a comma-separated key=value string. refers to the metadata in the asociated MANIFEST file.')
    cat('features: comma-separated list of features to plot along the x-axis. features are plotted in the order supplied.\n')
    stop("usage: consolidate_bed_enrichment_tables.R tag signal features out.csv enrich_object1.rda [ enrich_object2.rda ... enrich_objectN.rda ]")
}

group.tag <- args[1]
signal <- args[2]
features_string <- args[3]
outcsv <- args[4]
enrich.files <- args[-(1:4)]


if (file.exists(outcsv))
    stop(paste('output file', outcsv, 'already exists, please delete it first'))

suppressMessages(library(data.table))
suppressMessages(library(mutenrich))


signal.selections <- do.call(rbind, strsplit(strsplit(signal, ',')[[1]], '='))
colnames(signal.selections) <- c('key', 'value')
cat('selecting signal via\n')
print(signal.selections)

feature.selections <- strsplit(features_string, ',')[[1]]
cat('selecting features for plotting (in order):', feature.selections, '\n')

# Each .rda file contains:
#    - es    - either a single enrichment object or a list of enrichment objects.
#              If a list, names(es) are the feature names.
#    - emeta - data.frame with 1 row per feature and 1 column per metadata entry.
#              metadata entries are user defined, so no assumption should be made
#              about what columns are in this data.frame.

eslist <- lapply(enrich.files, function(efile) {
    loaded.objs <- load(efile)
    if (!('es' %in% loaded.objs & 'emeta' %in% loaded.objs))
        stop(paste('enrichment summary file', efile, 'expected to contain objects "es" and "emeta", but got', loaded.objs))

    if (!('list' %in% class(es))) {
        cat('ignoring signal: enrichment files contain a single signal\n')
        signal <- 'no signal selection'
        es[feature.selections,]  # features: reorder and subset
    } else {
        emeta.idxs <- rep(TRUE, nrow(emeta))
        for (i in 1:nrow(signal.selections)) {
            cat('selecting:', signal.selections[i,1], '=', signal.selections[i,2],'\n')
            cat('passing before:', sum(emeta.idxs))
            emeta.idxs <- emeta.idxs & (emeta[,signal.selections[i,1]] == signal.selections[i,2])
            cat(' passing after:', sum(emeta.idxs), '\n')
        }
        if (sum(emeta.idxs) > 1) {
            warning(paste('only one signal will be plotted but', sum(emeta.idxs), 'matched the signal criteria. using first matching signal'))
        }
        signal.to.keep <- which(emeta.idxs)[1]
        cat('keeping signal:', rownames(emeta)[signal.to.keep],'\n')
        es[[signal.to.keep]][feature.selections,]  # features: reorder and subset
    }
})

print(eslist)

# Write table of statistics
stat.table <- do.call(rbind, lapply(1:length(eslist), function(i) {
    es <- eslist[[i]]
    cbind(group.tag, class=rownames(es), es)
}))
fwrite(stat.table, file=outcsv)

if ('snakemake' %in% ls()) {
    sink()
}
