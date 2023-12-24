#!/usr/bin/env Rscript

# detect script being run by snakemake
# if so, make a mock commandArgs function
if ('snakemake' %in% ls()) {
    logfile <- snakemake@log[[1]]
    con <- file(logfile, 'w')
    sink(con, type='output')
    sink(con, type='message')

    commandArgs <- function(...) unlist(c(
        snakemake@input['metadata'],
        snakemake@output['mrca_timings'],
        snakemake@input[['mrcas']]   # list
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}


args <- commandArgs(trailingOnly=T)

if (length(args) < 3)
    stop("usage: finalize_mrca_timings.R metadata.csv out.timings.csv mrca1.txt [ mrca2.txt ... mrcaN.txt ]")

meta.file <- args[1]
out.timings <- args[2]
mrca.files <- args[-(1:2)]

for (f in c(out.timings)) {
    if (file.exists(f))
        stop(sprintf("output file %s already exists, please delete it first", f))
}

suppressMessages(library(data.table))

meta <- fread(meta.file)

# Ignore the SBS1 burden models in favor of the total mutation burden model
mrcas <- lapply(mrca.files, function(f) fread(f)[model == 'Total burden'])

z <- cbind(do.call(rbind, lapply(mrcas, function(x) {
    donor <- meta[sample == x$sample[1]]$donor[1]
    age <- meta[sample == x$sample[1]]$age[1]
    forward.estimate <- min(x$years.shared)
    min.private <- min(x$years.private)
    data.frame(donor=donor,
               pair=paste0(x$sample[1], ' vs. ', x$sample[2]),
               age=age,
               # forward/reverse method no longer used; kept for posterity
               forward.estimate=forward.estimate,
               min.private=min.private,
               reverse.estimate=pmax(forward.estimate, age - min.private),
               # these bounds are calculated using the 95% CIs on the linear model's
               # intercept and age slope.
               shared.est=mean(x$years.shared),
               shared.lb=min(x$years.shared.lb),
               shared.ub=max(x$years.shared.ub)
    )
})))
    
fwrite(z, file=out.timings)

if ('snakemake' %in% ls()) {
    sink()
}
