#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 4)
    stop('usage: gene_odds_ratio.R tumor.gene_mutations.txt neuron.gene_mutations.txt oligo.gene_mutations.txt out.rda')


tumor.txt <- args[1]
neuron.txt <- args[2]
oligo.txt <- args[3]
out.rda <- args[4]

if (file.exists(out.rda))
    stop(paste('output file', out.rda, 'already exists, please delete it first'))


data.cancer <- read.table(tumor.txt, sep='\t')
colnames(data.cancer) <- c('lineno', 'type', 'gene')
# when a single mut has multiple annotations, we want to keep all of the unique
# genes but keep at most one hit per gene.
data.cancer <- data.cancer[!duplicated(paste(data.cancer$lineno, data.cancer$gene)),]

data.sc <- c()
for (f in c(neuron.txt, oligo.txt)) {
    newdata <- read.table(f, sep='\t')
    colnames(newdata) <- c('lineno', 'type', 'gene')
    newdata <- newdata[!duplicated(paste(newdata$lineno, newdata$gene)),]
    data.sc <- rbind(data.sc, newdata)
}


types <- unique(data.cancer$type)
if (length(types) != 1)
    stop('this script only supports a single cancer type at a time')

top.max <- 10000
top.min  <- 1
tops <- seq(from=top.max,by=-10,to=top.min)

# there's only one type per run now
typelists <- lapply(types,function(type) {
  tmp <- table(data.cancer$gene[data.cancer$type == type])
  tmp <- names(tmp[order(tmp,decreasing=T)])
  return(tmp)
})[[1]]


cat('starting analysis\n')

topone <- function(top,tl,print=F,return.odds=F) {
  odds <- NA
  pval <- NA
  if(top <= length(tl)) {
    genes <- tl[1:top]
    tmp <- data.sc
    tmp$in.top <- tmp$gene %in% genes
    fishtab <- table(tmp[,c("type","in.top")])[c('neuron','oligo'),]
    if(print) {
        print(fishtab)
    }
    if (is.matrix(fishtab)) {
        tmp <- fisher.test(fishtab)
        if(is.infinite(tmp$estimate)) {
            tmp$estimate <- NA
        }
        odds <- tmp$estimate
        pval <- -log10(fisher.test(fishtab)$p.value)
    }
  }
  c(odds=odds, pval=pval)
}
  
cat('computing fisher tests\n')
fts <- sapply(tops, function(x) {
    topone(x, typelists)
})

cat("saving results to", out.rda, "\n")
save(types, tops, typelists, fts, file=out.rda)
