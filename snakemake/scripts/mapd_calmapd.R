cal_mapd <- function(ncnr){
	ncnr_1<-c(0,ncnr)
	n<-length(ncnr)
	b <- abs(log2(ncnr) - log2(ncnr_1[1:n]))
	return (median(b[2:n]))
}

dat <- read.table(f_in, header = T, sep="\t")
mapd <- cal_mapd(dat$cn.ratio)
cat (mapd)
