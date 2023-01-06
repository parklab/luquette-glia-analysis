#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)
if (!(length(args) %in% 4:5)) {
    cat('genes.bed is optional. if provided, only mutations in the gene symbols listed in genes.bed will be considered, and gene length info will be incorporated\n')
    stop('usage: gene_odds_ratio.R tumor.gene_mutations.txt neuron.gene_mutations.txt oligo.gene_mutations.txt genes.bed out.rda [ genes.bed ]')
}


tumor.txt <- args[1]
neuron.txt <- args[2]
oligo.txt <- args[3]
out.rda <- args[4]
genes.bed <- NULL
if (length(args) == 5)
    genes.bed <- args[5]

if (file.exists(out.rda))
    stop(paste('output file', out.rda, 'already exists, please delete it first'))

library(data.table)

# Handle the optional genes.bed. If not supplied, set to NULL here and make a dummy
# version based on the SnpEff input.
genes.info <- NULL
if (!is.null(genes.bed)) {
    genes.info <- fread(genes.bed)
    colnames(genes.info) <- c('chr', 'start', 'end', 'gene')
    genes.info[, length := end - start]
    setkey(genes.info, gene)
    # Get rid of duplicated gene names. The VAST majority of duplicated gene names
    # are "Y_RNA", "snoU13". There are 1751 duplicated gene names. Here they all are:
    #
    #> head(cumsum(sort(table(y[[4]]), decreasing=TRUE)),150)
    #          Y_RNA          snoU13              U3        SNORD112         SNORA70 
    #            747            1185            1257            1304            1331 
    #             U8         SNORA31         SNORA40         SNORA25         SNORA51 
    #           1357            1380            1399            1417            1430 
    #        SNORA63         SNORA48     Metazoa_SRP         SNORA26              U6 
    #           1443            1455            1466            1477            1488 
    #        SNORA75        SCARNA20         SNORA18         SNORA62         SNORA67 
    #           1496            1503            1510            1517            1524 
    #        SNORA72         SNORA73         SNORA74         snoU109          SNORA3 
    #           1531            1538            1545            1552            1558 
    #        SNORD56         SNORD74        SCARNA16        SCARNA18          SNORA1 
    #           1564            1570            1575            1580            1585 
    #        SNORA27         SNORA42         SNORA43         SNORA64         SNORA68 
    #           1590            1595            1600            1605            1610 
    #       SNORD113           ACA64        SCARNA11        SCARNA17        SCARNA21 
    #           1615            1619            1623            1627            1631 
    #        SNORA11         SNORA15         SNORA19         SNORA22         SNORA33 
    #           1635            1639            1643            1647            1651 
    #         SNORA4          SNORA7         SNORA77          SNORA8         SNORA81 
    #           1655            1659            1663            1667            1671 
    #        SNORD38         SNORD63         SNORD65         SNORD77         SNORD81 
    #           1675            1679            1683            1687            1691 
    #          snoZ6         5S_rRNA           ACA59         SNORA12          SNORA2 
    #           1695            1698            1701            1704            1707 
    #        SNORA20         SNORA30         SNORA32         SNORA57         SNORA76 
    #           1710            1713            1716            1719            1722 
    #         SNORA9        SNORD116         SNORD37         SNORD39         SNORD45 
    #           1725            1728            1731            1734            1737 
    #        SNORD46          SNORD5         SNORD60         SNORD78           Vault 
    #           1740            1743            1746            1749            1752 
    #          snR65 snoMe28S-Am2634         snoR442             7SK         ALG1L9P 
    #           1755            1758            1761            1763            1765 
    #        C2orf61        CYB561D2      DNAJC9-AS1           ELFN2         GOLGA8M 
    #           1767            1769            1771            1773            1775 
    #    IGHVII-31-1       LINC00484       LINC01115       LINC01238       LINC01347 
    #           1777            1779            1781            1783            1785 
    #      LINC01422       LINC01481       LINC01598           LYNX1            MAL2 
    #           1787            1789            1791            1793            1795 
    #      MIR3179-3         NBPF13P         OR7E47P       PROX1-AS1      RAET1E-AS1 
    #           1797            1799            1801            1803            1805 
    #           RGS5         RPS23P5        SCARNA15        SCARNA24         SCARNA4 
    #           1807            1809            1811            1813            1815 
    #        SCARNA6         SNORA17         SNORA21         SNORA24         SNORA36 
    #           1817            1819            1821            1823            1825 
    #        SNORA38         SNORA41         SNORA46         SNORA50         SNORA58 
    #           1827            1829            1831            1833            1835 
    #        SNORA66         SNORA71         SNORA79         SNORD11        SNORD111 
    #           1837            1839            1841            1843            1845 
    #       SNORD115         SNORD18         SNORD19        SNORD19B          SNORD2 
    #           1847            1849            1851            1853            1855 
    #        SNORD23         SNORD27         SNORD28         SNORD33         SNORD36 
    #           1857            1859            1861            1863            1865 
    #        SNORD41         SNORD42         SNORD43         SNORD51         SNORD59 
    #           1867            1869            1871            1873            1875 
    #        SNORD66         SNORD67         SNORD70         SNORD75         SPATA13 
    #           1877            1879            1881            1883            1885 
    #             U1              U4              U7     snoMBII-202        snoU2_19 
    #           1887            1889            1891            1893            1895 
    #   snosnR60_Z15            A1BG        A1BG-AS1            A1CF             A2M 
    #           1897            1898            1899            1900            1901 
    # A1BG and on are not duplicated (the count does not increase by >1).
    unique.genes <- setdiff(genes.info$gene, genes.info[duplicated(gene),gene])
    cat(sprintf('discarding %d out of %d genes for not being unique\n',
        nrow(genes.info)-length(unique.genes), nrow(genes.info)))
    genes.info <- genes.info[unique.genes]
    # Just combine non-unique gene names for now
}

data.cancer <- read.table(tumor.txt, sep='\t')
colnames(data.cancer) <- c('lineno', 'type', 'gene')
# when a single mut has multiple annotations, we want to keep all of the unique
# genes but keep at most one hit per gene.
data.cancer <- data.cancer[!duplicated(paste(data.cancer$lineno, data.cancer$gene)),]

data.sc <- c()
typenames <- c()
for (f in c(neuron.txt, oligo.txt)) {
    newdata <- read.table(f, sep='\t')
    colnames(newdata) <- c('lineno', 'type', 'gene')
    #newdata <- newdata[!duplicated(paste(newdata$lineno, newdata$gene)),]
    # new strategy: only take the first SnpEff hit. they are prioritized by functional
    # severity.  this is motivated by some gene clusters, particularly PCDHG{A,B}*,
    # that dominate the gene hit lists due to a large number of overlapping transcripts.
    newdata <- newdata[!duplicated(newdata$lineno),]
    # gene names in data.sc are XXX/YYY, where XXX usually = YYY
    newdata$gene <- sapply(strsplit(newdata$gene, '/'), head, 1)
    data.sc <- rbind(data.sc, newdata)
    typenames <- c(typenames, newdata$type[1])
}
cat('got typenames:\n')
print(typenames)


types <- unique(data.cancer$type)
if (length(types) != 1)
    stop('this script only supports a single cancer type at a time')

top.max <- 10000
top.min  <- 1
tops <- seq(from=top.max,by=-10,to=top.min)

# there's only one type per run now
typelists <- lapply(types,function(type) {
  # gene names in data.cancer are XXX/YYY, where XXX usually = YYY
  counts <- table(data.cancer$gene[data.cancer$type == type])
  # Faster to modify the names of the tabulated values rather than the data.cancer
  # data.frame since the data.frame has 1 row per mutation.
  genes <- sapply(strsplit(names(counts), '/'), head, 1)

  # no genes.bed supplied
  if (is.null(genes.info)) {
    genes.info <- data.table(gene=unique(genes), length=1)
    setkey(genes.info, gene)
  }

  names(counts) <- genes
  genes <- genes[genes %in% genes.info$gene]
  counts <- counts[genes]
  counts <- setNames(as.integer(counts), genes)  # getting rid of the annoying 'table' class
str(counts)
  genes.info[genes, count := counts[genes]]  # just joining the count for each gene
print(genes.info)
  genes.info[, norm.count := count / length]
print(genes.info)
  # XXX: normalizing by length was exploratory
  #return(genes.info[order(norm.count, decreasing=TRUE), gene])
  return(genes.info[order(count, decreasing=TRUE), gene])
})[[1]]


cat('starting analysis\n')

topone <- function(top,tl,print=F,return.odds=F) {
  odds <- NA
  pval <- NA
  if(top <= length(tl)) {
    genes <- tl[1:top]
    tmp <- data.sc
    tmp$in.top <- tmp$gene %in% genes
    fishtab <- table(tmp[,c("type","in.top")])[typenames,] # [c('neuron','oligo'),]
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
save(types, tops, typelists, fts, genes.info, file=out.rda)
