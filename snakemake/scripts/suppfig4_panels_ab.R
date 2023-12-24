#!/usr/bin/env Rscript

# detect script being run by snakemake
# if so, make a mock commandArgs function
if ('snakemake' %in% ls()) {
    logfile <- snakemake@log[[1]]
    con <- file(logfile, 'w')
    sink(con, type='output')
    sink(con, type='message')

    commandArgs <- function(...) unlist(c(
        snakemake@input['meta'],
        snakemake@input['matrix'],
        snakemake@output['panel_a_csv'],
        snakemake@output['panel_a_pdf'],
        snakemake@output['panel_b_csv'],
        snakemake@output['panel_b_pdf']
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}


args <- commandArgs(trailingOnly=T)

if (length(args) != 6)
    stop(sprintf("usage: suppfig4.R sample_metadata.csv SBS384.matrix.txt out.panel_a.csv out.panel_a.pdf out.panel_b.csv out.panel_b.pdf"))

meta.file <- args[1]
mat.file <- args[2]
panel.a.csv <- args[3]
panel.a.pdf <- args[4]
panel.b.csv <- args[5]
panel.b.pdf <- args[6]

for (f in c(panel.a.csv, panel.a.pdf, panel.b.csv, panel.b.pdf)) {
    if (file.exists(f))
        stop(sprintf("output file %s already exists, please delete it first", f))
}

# DO WORK
suppressMessages(library(data.table))
suppressMessages(library(scan2))

meta <- fread(meta.file)

x <- fread(mat.file)
x$tx <- substr(x[[1]], 1, 1)
x$mt <- substr(x[[1]], 5, 7)
x <- x[tx %in% c('T', 'U')]

y <- rbindlist(lapply(2:(ncol(x)-2), function(i) {
#y <- rbindlist(lapply(2, function(i) {
    this.idxs <- c(i, ncol(x), ncol(x)-1)
    ret <- cbind(sample=colnames(x)[i], x[,..this.idxs])
    colnames(ret)[2] <- 'count'
    # Need to sum over the 16 entries for each muttype, e.g. A[C>A]A, A[C>A]C, ...
    ret[,.(count=sum(count)),by=.(sample,mt,tx)]
}))

z <- meta[y,,on=.(sample)]


colors <- sbs96.cols[seq(1,96,16)]
names(colors) <- c('C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G')


# For Wilcoxon Rank Sum tests
op <- dcast(z[amp=='PTA' & type=='oligo'], amp + type + sample + age + mt ~ tx, value.var='count')
ops <- sapply(names(colors), function(this.mt) {
    w <- as.matrix(op[mt == this.mt, .(T,U)])
    wilcox.test(w[,1], w[,2], paired=FALSE)$p.value
})


np <- dcast(z[amp=='PTA' & type=='neuron'], amp + type + sample + age + mt ~ tx, value.var='count')
nps <- sapply(names(colors), function(this.mt) {
    w <- as.matrix(np[mt == this.mt, .(T,U)])
    wilcox.test(w[,1], w[,2], paired=FALSE)$p.value
})

# Barplots: get counts for all PTA neurons and PTA oligos
# Ordering needs to be interleaved so that T/U are displayed side by side
#    amp   type tx  mt    n
# 1: PTA neuron  T C>A  267
# 2: PTA neuron  T C>G  242
# 3: PTA neuron  T C>T  956
# 4: PTA neuron  T T>A  313
# 5: PTA neuron  T T>C 1796
# 6: PTA neuron  T T>G  304
# 7: PTA neuron  U C>A  221
# 8: PTA neuron  U C>G  465
# 9: PTA neuron  U C>T 1213
#10: PTA neuron  U T>A  349
#11: PTA neuron  U T>C  948
#12: PTA neuron  U T>G  439
nbp <- z[amp=='PTA' & type=='neuron',.(n=sum(count)),by=.(amp,type,tx,mt)][as.vector(t(matrix(1:12,ncol=2)))]
obp <- z[amp=='PTA' & type=='oligo',.(n=sum(count)),by=.(amp,type,tx,mt)][as.vector(t(matrix(1:12,ncol=2)))]

fwrite(rbind(nbp,obp), file=panel.a.csv)

pdf(file=panel.a.pdf, width=2.00, height=3, pointsize=5)
layout(1:2)
par(mar=c(4,2,2,1))
# returns x-locations of bars
bps <- barplot(nbp$n, names.arg=paste(nbp$tx, nbp$mt, sep=':'),
    col=colors[nbp$mt], border=NA, main='Neurons', las=3, ylim=c(0, max(nbp$n)*1.2))
bar.mids <- sapply(split(bps, rep(1:6, each=2)), mean)
# midpoint between each successive set of 2 bars
text(x=bar.mids, y=max(nbp$n)*1.1, labels=round(nps, 5), cex=0.8)

bps <- barplot(obp$n, names.arg=paste(nbp$tx, nbp$mt, sep=':'),
    col=colors[nbp$mt], border=NA, main='Oligodendrocytes', las=3, ylim=c(0, max(obp$n)*1.2))
bar.mids <- sapply(split(bps, rep(1:6, each=2)), mean)
text(x=bar.mids, y=max(obp$n)*1.1, labels=round(ops, 5), cex=0.8)
dev.off()

fwrite(z[,.(donor,sample,amp,type,age,mt,tx,count)], file=panel.b.csv)

pdf(width=6.5, height=3, pointsize=5, file=panel.b.pdf)
layout(matrix(1:12, ncol=6, byrow=T))
for (ct in c('neuron', 'oligo')) {
    for (this.mt in names(colors)) {
        zzt <- z[amp=='PTA' & type==ct & mt == this.mt & tx=='T']
        zzu <- z[amp=='PTA' & type==ct & mt == this.mt & tx=='U']
        par(mar=c(2,2,2,1))
        plot(c(zzt$age, zzu$age), c(zzt$count, zzu$count),
            pch=rep(c(1,16), each=nrow(zzt)), col=colors[this.mt],
            main=this.mt, xlab='Age', ylab='sSNVs')
        abline(lm(count~age, data=zzt), lty='solid', lwd=1)
        abline(lm(count~age, data=zzu), lty='dashed', lwd=1)
    }
}
dev.off()

if ('snakemake' %in% ls()) {
    sink()
}
