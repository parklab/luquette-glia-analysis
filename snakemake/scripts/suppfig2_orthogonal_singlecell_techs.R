#!/usr/bin/env Rscript

# detect script being run by snakemake
# if so, make a mock commandArgs function
if ('snakemake' %in% ls()) {
    logfile <- snakemake@log[[1]]
    con <- file(logfile, 'w')
    sink(con, type='output')
    sink(con, type='message')

    commandArgs <- function(...) unlist(c(
        snakemake@input['xlsx'],   # this contains NanoSeq and META-CS neuron results
        snakemake@output['pdf']
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}


args <- commandArgs(trailingOnly=T)

if (length(args) != 2) {
    cat("WARNING: this script uses an Excel table with pasted-in values for our PTA data and both external datasets. MAKE SURE TO UPDATE THIS MANUALLY!\n")
    stop("usage: suppfig2_orthogonal_singlecell_techs.R xlsx out.pdf")
}

duplex.xlsx <- args[1]
out.pdf <- args[2]

for (f in c(out.pdf)) {
    if (file.exists(f))
        stop(sprintf("output file %s already exists, please delete it first", f))
}


# DO WORK
library(readxl)

aging <- read_excel(duplex.xlsx, sheet='Aging rates')
snvs <- read_excel(duplex.xlsx, sheet='SNV counts per cell')
# Group by rough age (all within 1 year)
snvs$Group <- c(rep(1:3, c(10,11,11)), rep(1:3, each=3))
sigs <- read_excel(duplex.xlsx, sheet='SNV mutation signatures')


pdf(file=out.pdf, width=7, height=2.5, pointsize=8)

par(bty='l', mar=c(5,4,2.0,2))
layout(matrix(c(1,1,2,2, 3,3,3,3, 4,4,4,4, 5:7, 8), nrow=4), widths=c(1,1.75,1.5,3.5))

snv.rates <- as.numeric(aging[1,-1])
plot(snv.rates, pch=16, xlim=c(0.5,3.5),
    xaxt='n', xlab='', ylab='SNVs gained per year', ylim=c(12,20.5),
    main='SNV mutation rate')
axis(side=1, at=1:3, labels=colnames(aging[,-1]), las=3)
for (i in 2:4) {
    arrows(x0=i-1, y0=as.numeric(aging[2,i]), y1=as.numeric(aging[3,i]),
        length=0, lwd=1.5)
}

indel.rates <- as.numeric(aging[4,-1])
plot(indel.rates, pch=16, xlim=c(0.5,3.5),
    xaxt='n', xlab='', ylab='Indels gained per year', ylim=c(1,4),
    main='Indel mutation rate')
text(2.9, 2.5, 'N/A')
axis(side=1, at=1:3, labels=colnames(aging[,-1]), las=3)
for (i in 2:4) {
    arrows(x0=i-1, y0=as.numeric(aging[5,i]), y1=as.numeric(aging[6,i]),
        length=0, lwd=2)
}

par(bty='l', mar=c(5,4,3,2))

boxplot(SNVs.per.cell ~ Technology + Group, data=snvs, col=c('powderblue', 'orange'), ylim=c(250,1650), ylab='Somatic SNVs per single neuron', main='Total mutation burden per neuron', xaxt='n', xlab='', las=3, lty='solid')
stripchart(SNVs.per.cell ~ Technology + Group, data=snvs, vertical=TRUE, pch=20, method='jitter', add=TRUE)
axis(side=1, at=1:6, rep(c('META-CS', 'PTA'), 3), las=3)
arrows(x0=c(1,3,5)-0.2, x1=c(2,4,6)+0.2, y0=1600, code=0)
text(x=c(1,3,5)+0.5, y=1660, c('19 y.o.', '49 y.o.', '76-77 y.o.'))
legend('bottomright', fill=c('powderblue','orange'), legend=c('META-CS', 'PTA'))



# blank column
plot(1, pch=NA, bty='n', xaxt='n', yaxt='n', xlab='', ylab='')


# Something weird happened when pasting the data into Excel. Need to
# delete a trailing space, but no space character I know of matches it.
nanoseq <- as.numeric(substr(unlist(sigs[,2]),1,nchar(unlist(sigs[,2]))-1))
metacs <- as.numeric(substr(unlist(sigs[,3]),1,nchar(unlist(sigs[,3]))-1))
pta <- as.numeric(substr(unlist(sigs[,4]),1,nchar(unlist(sigs[,4]))-1))
mutsigs <- cbind(NanoSeq=nanoseq, `META-CS`=metacs, PTA=pta)

sbs96.cols <- rep(c("deepskyblue", "black", "firebrick2", "grey", "chartreuse3", "pink2"), each=16)

plotfun <- function(x, ylab, ...) {
    cossim <- function(a, b) round(sum(a*b) /(sqrt(sum(a^2))*sqrt(sum(b^2))), 3)
    p <- barplot(x/sum(x), space=0.5, border=NA, ylab='', xlab='', col=sbs96.cols, ...)
    abline(v = (p[seq(4, length(p) - 1, 4)] + p[seq(5, length(p),4)])/2, col = "grey")
    mtext(ylab, side=2, line=3, las=1, cex=4/5)

    # Don't plot the cossim=1 case
    sims <- c(cossim(x, nanoseq), cossim(x, metacs), cossim(x, pta))
    legend('topleft', title='Cosine sim.', bg='white', #bty='n',
        legend=paste(legend.text[sims != 1], sims[sims != 1]))
}

legend.text <- c('NanoSeq', 'META-CS', 'PTA')

par(mar=c(0,8,1,1))
plotfun(nanoseq, ylab='NanoSeq')
plotfun(metacs, ylab='META-CS')
plotfun(pta, ylab='PTA')

dev.off()


if ('snakemake' %in% ls()) {
    sink()
}
