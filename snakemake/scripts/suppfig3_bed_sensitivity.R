#!/usr/bin/env Rscript

# detect script being run by snakemake
# if so, make a mock commandArgs function
if ('snakemake' %in% ls()) {
    logfile <- snakemake@log[[1]]
    con <- file(logfile, 'w')
    sink(con, type='output')
    sink(con, type='message')

    commandArgs <- function(...) unlist(c(
        snakemake@output['pdf'],
        snakemake@params['datasources'],
        snakemake@params['group_colors'],
        snakemake@input[['sens']]
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}


args <- commandArgs(trailingOnly=T)

if (length(args) < 4) {
    cat("datasource_string is a comma-separated set of recognized datasources\n")
    cat("color_string is a comma-separated, key=value string mapping groups to colors\n")
    cat('e.g., pta_neuron=black,pta_oligo=red\n')
    stop("usage: suppfig3_bed_sensitivity.R out.pdf datasource_string color_string group1.summary.txt [ group2.summary.txt ... groupN.summary.txt ]")
}

outfile <- args[1]
datasource.string <- args[2]
color.string <- args[3]
sens.files <- args[-(1:3)]

for (f in c(outfile)) {
    if (file.exists(f))
        stop(sprintf("output file %s already exists, please delete it first", f))
}

# map group -> color
color.strings <- strsplit(color.string, split=',', fixed=TRUE)[[1]]
color.map <- setNames(sapply(strsplit(color.strings, split='=', fixed=TRUE), tail, 1),
                      sapply(strsplit(color.strings, split='=', fixed=TRUE), head, 1))
print(color.map)

# DO WORK
suppressMessages(library(scan2))

# gencode_simplified is in an earlier figure
#bed.dsrcs <- c('gencode_simplified', 'gtex', 'nott_enhprom', 'roadmap')
recognized.dsrcs <- c('roadmap', 'nott_enhprom', 'gtex')
requested.dsrcs <- strsplit(datasource.string, split=',', fixed=TRUE)[[1]]
if (!all(requested.dsrcs %in% recognized.dsrcs)) {
    cat('datasource_string: some requested datasources are not recognized:\n')
    print(requested.dsrcs)
    cat('recognized datasources (CASE SENSITIVE):\n')
    print(recognized.dsrcs)
    stop('see above error')
}
bed.dsrcs <- intersect(requested.dsrcs, recognized.dsrcs)

exclude.dcls <- c('genes')
exclude.features <- c('outside', '.', '8_ZNF/Rpts', '10_TssBiv', '11_BivFlnk', '12_EnhBiv', '3_TxFlnk')

# Unlike deciles that are easily ordered by numerical value, need to explicitly
# order the BED regions.
ordered.features <- c("1_TssA","2_TssAFlnk","4_Tx","5_TxWk","6_EnhG","7_Enh","8_ZNF/Rpts","9_Het","13_ReprPC","14_ReprPCWk","15_Quies","astrocyte_enhancer","neuron_enhancer","oligo_enhancer","microglia_enhancer","astrocyte_promoter","neuron_promoter","oligo_promoter","microglia_promoter","tx","utx")

# Exclude QBED records because there's a datasource=gtex QBED type that deals with expression
x <- rbindlist(lapply(sens.files, function(f) {
    ret <- fread(f)
    # MDA cells are not eligible for SCAN2 rescue
    if (substr(ret$group[1], 1, 4) == 'mda_') {
        this.passtype <- 'pass'
    } else {
        this.passtype <- 'rescue'
    }
    ret[bedtype == 'bed' & passtype == this.passtype & datasource %in% bed.dsrcs & !(dataclass %in% exclude.dcls) & !(feature %in% exclude.features)]
}))

x$feature <- factor(x$feature, levels=ordered.features, ordered=TRUE)

# data.table's integer64 causes as.matrix() to fail
x[, width := as.numeric(width)]


# for barplots, require.zero=TRUE should almost always be used
pf <- function(ys, col.map=NULL, require.zero=TRUE, do.layout=FALSE, show.xaxis=FALSE) {
    if (do.layout)
        layout(t(1:length(ys)))
    yvals <- as.matrix(do.call(cbind, lapply(ys, function(y) y[,-1])))
    if (require.zero)
        yvals <- c(0, yvals)
    ylim <- range(pretty(yvals))
    for (i in 1:length(ys)) {
        old.mar <- par(mar=c(1,3,0,0))
        y <- ys[[i]]
print(y)
        x <- y$feature  # region names
print(x)
        mat <- t(as.matrix(y[,-1]))
print(mat)

        if (is.null(col.map))
            col <- 'black'
        else
            col=col.map[rownames(mat)]
print(col)

        barplot(mat, names.arg=x, border=NA, beside=TRUE,
            col=col,
            ylim=ylim, bty='n', ylab='', xlab='', las=2,
            yaxt=ifelse(i==1, 's', 'n'), xaxt='n')
        par(mar=old.mar)
        if (show.xaxis) {
            old.mar <- par(mar=c(10,3,0,0))
            zero.mat <- matrix(0, ncol=ncol(mat), nrow=nrow(mat))
print(zero.mat)
            barplot(zero.mat, names.arg=x, beside=TRUE, border=NA,
                bty='n', ylab='', xlab='', las=2, yaxt='n')
            par(mar=old.mar)
        }
    }
}



n.groups <- length(sens.files)
n.features <- sapply(bed.dsrcs, function(ds) length(unique(x[datasource == ds]$feature)))

pdf(file=outfile, width=1/4 + sum(n.features)/12*n.groups, height=5, pointsize=5)

# last row is for x-axis labels with a dummy (blank) barplot to position them.
# needs to be interleaved, unlike other rows, for show.xaxis logic in pf().
# they will probably  need to be dragged up in an SVG editor later
lmat <- matrix(1:((6-1)*length(bed.dsrcs)), byrow=TRUE, ncol=length(bed.dsrcs))
lmat <- rbind(lmat, max(lmat) + matrix(1:(2*length(bed.dsrcs)), nrow=2))
print(lmat)

layout(lmat,
    # width of each column proportional to number of features
    widths=1.0+n.features,
    # give final row (x-axis) extra room for long labels
    heights=c(rep(1,6), 1.25))

# Decile genomic sizes (in MB)
ys <- setNames(lapply(bed.dsrcs, function(ds)
    dcast(x[datasource == ds, .(feature, dataclass, width=width/1e6)], feature ~ dataclass, value.var='width', fun.aggregate=function(x) head(x,1), fill=0, sep='|')
), bed.dsrcs)
pf(ys, require.zero=TRUE)

# _UN_weighted sequencing depth
ys <- setNames(lapply(bed.dsrcs, function(ds) {
    sc.dp <- dcast(x[datasource == ds], feature ~ group, value.var='unweighted.mean.sc.dp', fun.aggregate=function(x) head(x,1), fill=0, sep='|')
    sc.dp
}), bed.dsrcs)
pf(ys, col.map=color.map)


# _UN_weighted allele balance
ys <- setNames(lapply(bed.dsrcs, function(ds) dcast(x[datasource == ds], feature ~ group, value.var='unweighted.mean.abs.gp.mu', fun.aggregate=function(x) head(x,1), fill=0, sep='|')), bed.dsrcs)
pf(ys, col.map=color.map)

ys <- setNames(lapply(bed.dsrcs, function(ds) dcast(x[datasource == ds], feature ~ group, value.var='unweighted.mean.gp.sd', fun.aggregate=function(x) head(x,1), fill=0, sep='|')), bed.dsrcs)
pf(ys, col.map=color.map)

# Weighted SNV sensitivity
ys <- setNames(lapply(bed.dsrcs, function(ds) dcast(x[datasource == ds & muttype == 'snv'], feature ~ group, value.var='weighted.sens', sep='|')), bed.dsrcs)
pf(ys, col.map=color.map)

# Weighted indel sensitivity
ys <- setNames(lapply(bed.dsrcs, function(ds) dcast(x[datasource == ds & muttype == 'indel'], feature ~ group, value.var='weighted.sens', fun.aggregate=function(x) head(x,1), fill=0, sep='|')), bed.dsrcs)
pf(ys, col.map=color.map, show.xaxis=TRUE)

dev.off()

if ('snakemake' %in% ls()) {
    sink()
}
