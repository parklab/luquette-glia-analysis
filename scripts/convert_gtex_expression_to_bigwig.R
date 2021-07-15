#!/usr/bin/env Rscript
#SBATCH -p park
#SBATCH -A park_contrib
#SBATCH -t 1:00:00
#SBATCH --mem=16G

args <- commandArgs(trailingOnly=TRUE)

if (length(args) != 6) {
    cat('the out_fmt.bed parameter requires a %s placeholder that will be\n')
    cat('replaced with quartile, quintile and decile.\n')
    stop('usage: make_gb.r gencode_gene_model.gtf tissue_with_spaces tissue_median_tpm.gct min.tpm out.rda out_fmt.bed')
}

gene.model.file <- args[1]
tissue <- args[2]
tissue.median.tpm.file <- args[3]
min.tpm <- as.numeric(args[4])
outrda <- args[5]
outfmtbed <- args[6]

if (length(grep('%s', outfmtbed)) == 0)
    stop('out_fmt.bed must include a %s placeholder')

if (file.exists(outrda))
    stop(paste('output file', outrda, 'already exists, please delete it first'))

for (ile in c('quartile', 'quintile', 'decile')) {
    outf <- sprintf(outfmtbed, ile)
    if (file.exists(outf))
        stop(paste('output file', outf, 'already exists, please delete it first'))
}


suppressMessages(library(GenomicRanges))
suppressMessages(library(rtracklayer))
suppressMessages(library(data.table))
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg19))
suppressMessages(library(bedtoolsr))
cat('Loading gene model:', gene.model.file, '\n')
# this version is the result of running GTEx's collapse script
gtex <- import.gff(gene.model.file)
# ENSEMBL IDs have underscores after them in GENCODE but not in GTEx
gtex$gene_id2 <- sapply(strsplit(gtex$gene_id,'_'),head,1)

cat('Annotating tissue median TPM expression levels:', tissue.median.tpm.file, '\n')
gtex.tpm <- fread(tissue.median.tpm.file, skip=2)
setkey(gtex.tpm, Name)
gtex.tpm <- gtex.tpm[gtex$gene_id2]
for (cn in colnames(gtex.tpm)[-(1:2)])
    mcols(gtex)[cn] <- gtex.tpm[[cn]]

gtex$mean.expr <- rowMeans(gtex.tpm[,-(1:2)])
e <- gtex.tpm[[tissue]]

# Only retain genes that have >1 TPM on average in either all samples
# or just brain specific samples.
# This helps to differentiate between low expression genes that are
# likely real genes and low expression regions of the genome that
# are more likely to be noise or spurious expression.
selected.genes <- as.character(seqnames(gtex)) %in% paste0('chr', 1:22) &
    (!is.na(gtex$mean.expr) & (gtex$mean.expr >= min.tpm | e >= min.tpm)) &
    apply(gtex.tpm[,-(1:2)],1,max) > 0
    #gtex$gene_type == 'protein_coding'

cat(sprintf('    %d / %d genes in autosomes, non-NA expression and median expression > 0 in at least 1 of %d tissues\n',
    sum(selected.genes), length(selected.genes), ncol(gtex.tpm)-2))

gtex <- gtex[selected.genes,]
e <- e[selected.genes]

# Set hg19 chr1-22 lengths for tileGenome
seqlevels(gtex) <- seqlevels(gtex)[1:22]
seqinfo(gtex) <- seqinfo(BSgenome.Hsapiens.UCSC.hg19)[paste0('chr',1:22),]
cat('Tiling hg19 with 100 bp bins\n')
bins <- tileGenome(seqlengths=seqlengths(gtex), tilewidth=100, cut.last.tile.in.chrom=T)

cat('Mapping GTEx to non-overlapping 100 bp bins and taking max expression\n')
ols <- findOverlaps(bins, gtex)
system.time(maxexpr <- sapply(split(e[to(ols)], from(ols)), max))
emap <- data.frame(from=as.integer(names(maxexpr)), to=maxexpr)
bins$max.expr <- NA
bins$max.expr[emap$from] <- emap$to

# Map expression levels to various quantiles.
cat('Mapping expression levels to various quantiles\n')
# quartiles
bins$qq <- Rle(findInterval(bins$max.expr,
    quantile(bins$max.expr, na.rm=T, probs=1:4/4),
    rightmost.closed=T)+1)
bins$qrunid <- Rle(rep(1:nrun(bins$qq), runLength(bins$qq)))
print(table(bins$qq))
# quintiles
bins$pq <- Rle(findInterval(bins$max.expr,
    quantile(bins$max.expr, na.rm=T, probs=1:5/5),
    rightmost.closed=T)+1)
bins$prunid <- Rle(rep(1:nrun(bins$pq), runLength(bins$pq)))
print(table(bins$pq))
# deciles
bins$dq <- Rle(findInterval(bins$max.expr,
    quantile(bins$max.expr, na.rm=T, probs=1:10/10),
    rightmost.closed=T)+1)
bins$drunid <- Rle(rep(1:nrun(bins$dq), runLength(bins$dq)))
print(table(bins$dq))
bins$lq <- Rle(findInterval(bins$max.expr,
    quantile(bins$max.expr, na.rm=T, probs=1:50/50),
    rightmost.closed=T)+1)
bins$lrunid <- Rle(rep(1:nrun(bins$lq), runLength(bins$lq)))
print(table(bins$lq))

print(bins)

# dataframe so we can use write.table
d <- as.data.frame(bins[!is.na(bins$qq),])
d$seqnames <- as.character(d$seqnames)
print(head(d))
for (ile in c('quartile', 'quintile', 'decile', '50ile')) {
    outf <- sprintf(outfmtbed, ile)
    select <- c(quartile='q', quintile='p', decile='d', `50ile`='l')[ile]
    
    tmpbed <- tempfile(tmpdir='/n/data1/hms/dbmi/park/jluquette/pta/gtex/tmp', fileext='.bed')
    cat(paste("Writing ungrouped BED to", tmpbed, "\n"))
    write.table(d[,c('seqnames','start','end',paste0(select,'q'),paste0(select,'runid'))],
        file=tmpbed, quote=F, sep='\t', row.names=F, col.names=F)

    cat("Grouping bed by run\n")
    # first two columns re-iterate the group by criteria
    # means are unnecessary, just wanted to see if any non-integers
    # were returned, which would indicate a bug.  all runs should have
    # identical q and runid values.
    dg <- bt.groupby(i=tmpbed, g='1,5', c='1,2,3,4,5',
        o='first,first,last,mean,mean')[,-(1:2)]
    dg[,1] <- as.character(dg[,1])
    print(head(dg))

    cat("Writing final grouped BED\n")
    write.table(dg, file=outf, quote=F, row.names=F, col.names=F, sep='\t')
}

save(gtex, e, bins, file=outrda, compress=FALSE)

cat("Done\n")
print(gc())
