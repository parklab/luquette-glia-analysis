#!/usr/bin/env Rscript

# detect script being run by snakemake
# if so, make a mock commandArgs function
if ('snakemake' %in% ls()) {
    logfile <- snakemake@log[[1]]
    con <- file(logfile, 'w')
    sink(con, type='output')
    sink(con, type='message')

    commandArgs <- function(...) unlist(c(
        snakemake@input['muts'],
        snakemake@input['meta'],
        snakemake@output['csv']
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}


args <- commandArgs(trailingOnly=T)

if (length(args) != 3)
    stop("usage: supptable2.R allmuts_unfiltered.csv metadata.csv out.csv")

muts.file <- args[1]
meta.file <- args[2]
outfile <- args[3]

for (f in c(outfile)) {
    if (file.exists(f))
        stop(sprintf("output file %s already exists, please delete it first", f))
}

# DO WORK
library(data.table)

muts <- fread(muts.file)
meta <- fread(meta.file)


muts <- meta[,.(sample,amp,type,batch)][muts,,on=.(sample)]

# Printing out this table just for diagnostics.
cat("Including unsupported SCAN2 calls:\n")
addmargins(table(muts[final.filter == FALSE, .(paste(muttype,rescue), paste(amp,type,batch))]))


# For PTA samples, all calls are used.
# For MDA samples, indel calls are not supported (though nothing stops SCAN2 from
# outputting them) and rescue is not applicable (though, again, SCAN2 signature
# rescue was run to create uniform output formats. The signature rescued calls
# were then discarded).
muts <- muts[amp != 'MDA' | (rescue == FALSE & muttype=='snv')]
cat("Excluding unsupported SCAN2 calls:\n")
addmargins(table(muts[final.filter == FALSE, .(paste(muttype,rescue), paste(amp,type,batch))]))

# Get rid of added columns after filtration is complete
muts$amp <- NULL
muts$type <- NULL
muts$batch <- NULL

fwrite(muts, file=outfile)

if ('snakemake' %in% ls()) {
    sink()
}
