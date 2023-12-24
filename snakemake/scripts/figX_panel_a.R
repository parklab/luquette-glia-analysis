#!/usr/bin/env Rscript

# detect script being run by snakemake
# if so, make a mock commandArgs function
if ('snakemake' %in% ls()) {
    logfile <- snakemake@log[[1]]
    con <- file(logfile, 'w')
    sink(con, type='output')
    sink(con, type='message')

    commandArgs <- function(...) unlist(c(
        snakemake@input['csv'],
        snakemake@output['pdf']
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}


args <- commandArgs(trailingOnly=T)

if (length(args) != 2) {
    stop(sprintf("usage: figX_panel_a.R mrca_timings.csv out.pdf"))
}

infile <- args[1]
outfile <- args[2]

for (f in c(outfile)) {
    if (file.exists(f))
        stop(sprintf("output file %s already exists, please delete it first", f))
}

# DO WORK
library(data.table)

mrca.data <- fread(infile)

par(mar=c(4,2,0.1,0.1))
#left.limit <- -30
left.limit <- -82
rect.size <- 1/12
lwd <- 4

# Row 1: left side of box (count-forward time to MRCA, may be negative, meaning pre-birth (birth=0))
# Row 2: right side of box (count-backward time to MRCA)
# Row 3: age of individual (endpoint of arrow)
# Column 1: MDA eldest pair: GliaLC-4-{F11,G10}
# Column 2: PTA eldest pair: 5657_OL_{4,6}
# Column 3: youngest pair: 5559-Oligo-{5,8}
# Note: rows are plotted from bottom-up
z <- rbind(
    #mrca.data$forward.estimate,
    #mrca.data$reverse.estimate,
    mrca.data$shared.lb,
    mrca.data$shared.ub,
    mrca.data$age,
    mrca.data$shared.est
)

pdf(width=6.5, height=2.5, pointsize=8, file=outfile)

development.col <- 'orange'
aging.col <- 'dodgerblue3'

plot(0, pch=NA, xlim=c(left.limit,82), ylim=c(0.5,3.5), bty='l', xlab='Time (years)', yaxt='n', xaxt='n')

for (i in 1:ncol(z)) {
    rect(xleft=z[1,i], xright=z[2,i], ybottom=i-rect.size, ytop=i+rect.size)
    arrows(y0=i,y1=i+c(1,-1)*rect.size*2,x0=z[2,i],x1=z[3,i],lwd=lwd,code=0, col=aging.col)
    segments(y0=i,y1=i,x0=left.limit+2,x1=min(0,z[1,i]), lwd=lwd, col=development.col)
    segments(y0=i,y1=i,x0=min(0,z[1,i]), x1=z[1,i], lwd=lwd, col=aging.col)
    points(x=left.limit+2,y=i,pch=20)
    points(x=z[4,i], y=i, pch=20)
}
abline(v=0, lty='dashed', col='darkgrey')
axis(side=1, at=c(left.limit+2, seq(0,80,20)), labels=c('Zygote', 'Birth', seq(20,80,20)))

dev.off()

if ('snakemake' %in% ls()) {
    sink()
}
