library("DNAcopy")

lowess.gc <- function(jtkx, jtky) {
        jtklow <- lowess(jtkx, log(jtky), f=0.05)
        jtkz <- approx(jtklow$x, jtklow$y, jtkx)
        return(exp(log(jtky) - jtkz$y))
}

cbs.segment01 <- function(f_in, f_out, f_out_short, bad.bins, varbin.gc, alpha, nperm, undo.SD, min.width) {
    #####################################
    # f_in: input file (output from varbin.50k.bam.py)
    # f_out: output file
    # bad.bins: path to bad bins file
    # varbin.gc: path to gc content file
    # sample.name: sample name
    ####################################

	gc <- read.table(varbin.gc, header=T)
	bad <- read.table(bad.bins, header=F)

	chrom.numeric <- substring(gc$bin.chrom, 4)
	chrom.numeric[which(gc$bin.chrom == "chrX")] <- "23"
	chrom.numeric[which(gc$bin.chrom == "chrY")] <- "24"
	chrom.numeric <- as.numeric(chrom.numeric)

	thisRatio <- read.table(f_in, header=F) 
	# thisRatio <- read.table(paste(indir, varbin.data, sep="/"), header=F) 
	names(thisRatio) <- c("chrom", "chrompos", "abspos", "bincount", "ratio")
	thisRatio$chrom <- chrom.numeric
	a <- thisRatio$bincount + 1
	#thisRatio$ratio <- a / mean(a)
	thisRatio$ratio <- a / median(a)
	thisRatio$gc.content <- gc$gc.content
	thisRatio$lowratio <- lowess.gc(thisRatio$gc.content, thisRatio$ratio)

	a <- quantile(gc$bin.length, 0.985)
	thisRatioNobad <- thisRatio[which(bad[, 1] == 0), ]
	
	set.seed(25) 
	CNA.object <- CNA(log(thisRatio$lowratio, base=2), thisRatio$chrom, thisRatio$chrompos, data.type="logratio", sampleid=tail(strsplit(f_out, '/'), n=1)) 
	smoothed.CNA.object <- smooth.CNA(CNA.object) 
	segment.smoothed.CNA.object <- segment(smoothed.CNA.object, alpha=alpha, nperm=nperm, undo.splits="sdundo", undo.SD=undo.SD, min.width=min.width) 
	thisShort <- segment.smoothed.CNA.object[[2]]

	m <- matrix(data=0, nrow=nrow(thisRatio), ncol=1)	
	prevEnd <- 0
	for (i in 1:nrow(thisShort)) {
		thisStart <- prevEnd + 1
		thisEnd <- prevEnd + thisShort$num.mark[i]
		m[thisStart:thisEnd, 1] <- 2^thisShort$seg.mean[i]
		prevEnd = thisEnd
	}
	
	thisRatio$seg.mean.LOWESS <- m[, 1]

	chr <- thisRatio$chrom
	chr.shift <- c(chr[-1], chr[length(chr)])
	vlines <- c(1, thisRatio$abspos[which(chr != chr.shift) + 1], thisRatio$abspos[nrow(thisRatio)])
	hlines <- c(0.5, 1.0, 1.5, 2.0)
	chr.text <- c(1:22, "X", "Y")
	chr.text <- chr.text[c(1:14, 16, 18, 20, 22:24)]
	vlines.shift <- c(vlines[-1], 4*10^9)
	chr.at <- vlines + (vlines.shift - vlines) / 2
	chr.at <- chr.at[c(1:14, 16, 18, 20, 22:24)]
	x.at <- c(0, 0.5, 1, 1.5, 2, 2.5, 3) * 10^9
	x.labels <- c("0", "0.5", "1.0", "1.5", "2.0", "2.5", "3.0")
	y.at <- c(0.005, 0.020, 0.100, 0.500, 2.000)
	y.labels <- c("0.005", "0.020", "0.100", "0.500", "2.000")

	write.table(thisRatio, sep="\t", file=f_out, quote=F, row.names=F) 
	write.table(thisShort, sep="\t", file=f_out_short, quote=F, row.names=F) 

}

args = commandArgs(TRUE)
file_in = args[1]
file_out = args[2]
file_out_short = args[3]
bad.bins = args[4]
varbin.gc = args[5]

cbs.segment01(file_in, file_out, file_out_short, bad.bins=bad.bins, varbin.gc=varbin.gc, alpha=0.02, nperm=1000, undo.SD=1.0, min.width=5)
