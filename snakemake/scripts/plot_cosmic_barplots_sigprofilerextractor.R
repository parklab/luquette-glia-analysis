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
        snakemake@input['expomat'],
        snakemake@input['metadata'],
        snakemake@output['svg'],
        snakemake@output['pdf'],
        snakemake@output['csv']
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 5) {
    stop("usage: plot_cosmic_barplots_sigprofilerextractor.R expomat.csv metadata.csv out.svg out.pdf out.csv")
}

expomat.csv <- args[1]
metadata.csv <- args[2]
out.svg <- args[3]
out.pdf <- args[4]
out.csv <- args[5]

if (file.exists(out.svg))
    stop(paste('output file', out.svg, 'already exists, please delete it first'))
if (file.exists(out.pdf))
    stop(paste('output file', out.pdf, 'already exists, please delete it first'))

suppressMessages(library(data.table))
suppressMessages(library(svglite))


E <- fread(expomat.csv)
head(E)

meta <- fread(metadata.csv)[sample %in% E$Samples]  # needed for age/sorting
meta <- meta[order(age)]

# Reorder E by increasing age
signames <- colnames(E)[-1]  # first column is Samples
setkey(E, Samples)
E <- E[meta$sample]   # Sort to match meta's order (sorted by age, above)
cat("Age sorted:\n")
head(E)
Emat <- as.matrix(E[,-1])
rownames(Emat) <- E$Samples


sigcols <-c(SBS5="#1f78b4", SBS1="#e31a1c", SBS32="#fdbf6f", SBS16="#33a02c", SBS19="grey")

# Sort the signatures from low to high contribution so that
# the largest contributors are the bars closest to the bottom
# of the plot.
# rearrange sigcols too so the same colors are used across plots
#sigcols <- sigcols[order(colSums(Emat), decreasing=T)]
signames <- signames[order(colSums(Emat), decreasing=T)]
Emat <- Emat[,order(colSums(Emat), decreasing=T)]
cat("Final reordered matrix:\n")
head(Emat)
cat("Column sums:\n")
print(colSums(Emat))

# get height (in lines) of longest sample name (will be in x-axis label)
em <- strheight('M', unit='inches')  # height of a capital M in inches

xheight <- max(strwidth(E$Samples, unit='inches'))/em  # size in lines, for par(mar)
# give 3.5 inches for the plot area, xheight for x-axis and 2 lines for title
# xheight*em = size in INCHES instead of lines
figheight <- 2*(3.5 + xheight*em + 2*em)

# 4 - y axis, 1 - right margin
figwidth <- (4 + 1)*em + nrow(E)*em*1.2  # in inches

fwrite(E, file=out.csv)

# Loop over devices to save both pdf and svgs
devs <- list(svglite, pdf)
outs <- c(out.svg, out.pdf)
for (i in 1:2) {
    devs[[i]](file=outs[i], width=figwidth, height=figheight) #, pointsize=5)  # pointsize works, but haven't settled yet on figure sizes
    layout(1:2)
    par(mar=c(xheight,4,1,1))
    barplot(t(Emat), col=sigcols[colnames(Emat),], las=3, border=NA, ylab='Number of mutations', cex.names=3/4)
    legend('topleft', legend=signames, fill=sigcols[colnames(Emat),])

    x <- apply(Emat, 1, function(col) col/sum(col))
    barplot(x, col=sigcols[colnames(Emat),], las=3, border=NA, ylab='Fraction of mutations', cex.names=3/4)
    dev.off()
}

if ('snakemake' %in% ls()) {
    sink()
}
