#!/usr/bin/env Rscript

# detect script being run by snakemake
# if so, make a mock commandArgs function
if ('snakemake' %in% ls()) {
    logfile <- snakemake@log[[1]]
    con <- file(logfile, 'w')
    sink(con, type='output')
    sink(con, type='message')

    commandArgs <- function(...) unlist(c(
        snakemake@input['neuron_sens'],
        snakemake@input['oligo_sens'],
        snakemake@output['pdf']
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}


args <- commandArgs(trailingOnly=T)

if (length(args) != 3)
    stop("usage: suppfig3_quantile_sensitivity.2.R neuron.summary.txt oligo.summary.txt out.pdf")

ns.file <- args[1]
os.file <- args[2]
outfile <- args[3]

for (f in c(outfile)) {
    if (file.exists(f))
        stop(sprintf("output file %s already exists, please delete it first", f))
}

# DO WORK
suppressMessages(library(scan2))

ns <- fread(ns.file)
os <- fread(os.file)

q.dsrcs <- c('scrnaseq', 'scatacseq', 'repliseq', 'active_histone_mark', 'inactive_histone_mark', 'gtex', 'cancer_snvdens')
# Myeloid-MDS - extremely low number of mutations (~2k) across all tumors. Created very
# few deciles as a result. Should probably remove from all analyses, not just this.
exclude.dcls <- c('Neurons', 'Endothelial', 'endothelial', 'genes', 'tx_or_utx', 'Myeloid-MDS')  # genes, tx_or_utx are gtex BED region datasources

x <- rbind(ns,os)[passtype == 'rescue' & datasource %in% q.dsrcs & !(dataclass %in% exclude.dcls) & (binsize==1000 | (binsize==1000000 & datasource=='cancer_snvdens')) & n.quantiles==10 & feature %in% 1:10]

# data.table's integer64 causes as.matrix() to fail
x[, width := as.numeric(width)]

# Reverse RepliSeq decile order so the low=early and high=late
x[datasource == 'repliseq', feature := as.character(11 - as.integer(feature))]

#light.cols <- c('#a6cee3','#b2df8a','#fb9a99','#fdbf6f','#cab2d6','#ffff99')
#dark.cols  <- c('#1f78b4','#33a02c','#e31a1c','#ff7f00','#6a3d9a','#b15928')

lty.map <- c(
`pta_neuron`='solid',
`pta_oligo`='dashed'
)

color.map <- c(
`Astrocytes`='#F4C013',
`Endothelial`='#65B648',
`Excitatory-Neurons`='#000000',
`Inhibitory-Neurons`='#52C0CC',
`Microglia`='#A8388C',
`OPCs`='#929393',
`Oligodendrocytes`='#E91D21',
`astrocyte`='#F4C013',
`endothelial`='#65B648',
`excitatory_neuron`='#000000',
`inhibitory_neuron`='#52C0CC',
`microglia`='#A8388C',
`OPC`='#929393',
`oligo`='#E91D21',
`BG02ES`="grey",
`BJ`="grey",
`GM06990`="grey",
`GM12801`="grey",
`GM12812`="grey",
`GM12813`="grey",
`GM12878`="grey",
`HUVEC`="grey",
`HeLa-S3`="grey",
`HepG2`="grey",
`IMR90`="grey",
`K562`="grey",
`MCF-7`="grey",
`NHEK`="grey",
`SK-N-SH`="grey",
`H3K27ac`='#66c2a5',
`H3K36me3`='#fc8d62',
`H3K4me1`='#8da0cb',
`H3K4me3`='#e78ac3',
`H3K9ac`='#a6d854',
`H3K27me3`='#000000',
`H3K9me3`='#e5c494',
`Adipose_-_Subcutaneous`="grey",
`Adipose_-_Visceral__Omentum_`="grey",
`Adrenal_Gland`="grey",
`Artery_-_Aorta`="grey",
`Artery_-_Coronary`="grey",
`Artery_-_Tibial`="grey",
`Bladder`="grey",
`Brain_-_Amygdala`="green",
`Brain_-_Anterior_cingulate_cortex__BA24_`="green",
`Brain_-_Caudate__basal_ganglia_`="green",
`Brain_-_Cerebellar_Hemisphere`="green",
`Brain_-_Cerebellum`="green",
`Brain_-_Cortex`="green",
`Brain_-_Frontal_Cortex__BA9_`="green",
`Brain_-_Hippocampus`="green",
`Brain_-_Hypothalamus`="green",
`Brain_-_Nucleus_accumbens__basal_ganglia_`="green",
`Brain_-_Putamen__basal_ganglia_`="green",
`Brain_-_Spinal_cord__cervical_c-1_`="green",
`Brain_-_Substantia_nigra`="green",
`Breast_-_Mammary_Tissue`="grey",
`Cells_-_Cultured_fibroblasts`="grey",
`Cells_-_EBV-transformed_lymphocytes`="grey",
`Cervix_-_Ectocervix`="grey",
`Cervix_-_Endocervix`="grey",
`Colon_-_Sigmoid`="grey",
`Colon_-_Transverse`="grey",
`Esophagus_-_Gastroesophageal_Junction`="grey",
`Esophagus_-_Mucosa`="grey",
`Esophagus_-_Muscularis`="grey",
`Fallopian_Tube`="grey",
`Heart_-_Atrial_Appendage`="grey",
`Heart_-_Left_Ventricle`="grey",
`Kidney_-_Cortex`="grey",
`Kidney_-_Medulla`="grey",
`Liver`="grey",
`Lung`="grey",
`Minor_Salivary_Gland`="grey",
`Muscle_-_Skeletal`="grey",
`Nerve_-_Tibial`="grey",
`Ovary`="grey",
`Pancreas`="grey",
`Pituitary`="grey",
`Prostate`="grey",
`Skin_-_Not_Sun_Exposed__Suprapubic_`="grey",
`Skin_-_Sun_Exposed__Lower_leg_`="grey",
`Small_Intestine_-_Terminal_Ileum`="grey",
`Spleen`="grey",
`Stomach`="grey",
`Testis`="grey",
`Thyroid`="grey",
`Uterus`="grey",
`Vagina`="grey",
`Whole_Blood`="grey",
`Biliary-AdenoCA`="grey",
`Bladder-TCC`="grey",
`Bone-Cart`="grey",
`Bone-Epith`="grey",
`Bone-Leiomyo`="grey",
`Bone-Osteosarc`="grey",
`Breast-AdenoCa`="grey",
`Breast-DCIS`="grey",
`Breast-LobularCa`="grey",
`CNS-GBM`="orange",
`CNS-Medullo`="black",
`CNS-Oligo`="red",
`CNS-PiloAstro`="purple",
`Cervix-AdenoCA`="grey",
`Cervix-SCC`="grey",
`ColoRect-AdenoCA`="grey",
`Eso-AdenoCa`="grey",
`Head-SCC`="grey",
`Kidney-ChRCC`="grey",
`Kidney-RCC`="grey",
`Liver-HCC`="grey",
`Lung-AdenoCA`="grey",
`Lung-SCC`="grey",
`Lymph-BNHL`="grey",
`Lymph-CLL`="grey",
`Lymph-NOS`="grey",
`Myeloid-AML`="grey",
`Myeloid-MDS`="grey",
`Myeloid-MPN`="grey",
`Ovary-AdenoCA`="grey",
`Panc-AdenoCA`="grey",
`Panc-Endocrine`="grey",
`Prost-AdenoCA`="grey",
`Skin-Melanoma`="grey",
`Stomach-AdenoCA`="grey",
`Thy-AdenoCA`="grey",
`Uterus-AdenoCA`="grey"
)


pf <- function(ys, col.map=1:6, lty.map='solid', require.zero=FALSE, do.layout=FALSE, show.xaxis=FALSE) {
    if (do.layout)
        layout(t(1:length(ys)))
    yvals <- as.matrix(do.call(cbind, lapply(ys, function(y) y[,-1])))
    if (require.zero)
        yvals <- c(0, yvals)
    ylim <- range(pretty(yvals))
    for (i in 1:length(ys)) {
        old.mar <- par(mar=c(2,3,0,0.2))
        y <- ys[[i]]
        x <- y[,1]
        mat <- as.matrix(y[,-1])
        # first value is always dataclass
        col <- col.map[sapply(strsplit(colnames(mat), '|', fixed=TRUE), head, 1)]
print(sapply(strsplit(colnames(mat), '|', fixed=TRUE), head, 1))
print(col)
        if (length(lty.map) == 1)
            lty <- lty.map
        else
            lty <- sapply(strsplit(colnames(mat), '|', fixed=TRUE), function(els) if (length(els)>1) lty.map[els[2]] else 'solid')
print(lty)
        # For RepliSeq, we show the average of 16 cell lines. There's
        # no interesting difference between them.
        if (names(ys)[i] == 'repliseq') {
            mat <- cbind(rowMeans(mat[, lty=='solid']), rowMeans(mat[, lty=='dashed']))
            col <- 'grey'
            lty <- c('solid', 'dashed')
        }

print(colnames(mat))
print(lty)
        # Use a thinner line for GTEx expression since there are 54 lines to plot per group
        # or cancer_snvdens with 37 lines.
        matplot(x, mat, type='l', lty=lty, col=col, lwd=ifelse(names(ys)[i] %in% c('gtex', 'cancer_snvdens'), 1/4, 2/3), las=1,
            #main=names(ys)[i], 
            ylim=ylim, bty='n', ylab='', xlab='',
            yaxt=ifelse(i==1, 's', 'n'), xaxt='n')
        # for gtex and cancer_snvdens only, replot the non-grey colored lines on top so they're visible
        if (names(ys)[i] %in% c('gtex', 'cancer_snvdens')) {
cat('------------ values -------------\n')
print(x)
print(mat)
            nongrey <- mat[,which(col!='grey')]
            cat('overplotting', ncol(nongrey), 'emphasized lines\n')
            matplot(x, nongrey, type='l', lty=lty[col!='grey'], col=col[col != 'grey'], lwd=1/4,
                bty='n', ylab='', xlab='', add=TRUE)
        }
        if (show.xaxis) {
            # For some reason, R refuses to plot the tick labels at 1,5,10 below.
            # Just let it do its default (which unfortunately starts at 2).
            axis(side=1)#, at=1:10, labels=c(1,'','','',5,'','','','',10)) 
        }
        par(mar=old.mar)
    }
}


pdf(file=outfile, width=5.5, height=5, pointsize=6)
# last row is for color legends
layout(matrix(1:((6+1)*length(q.dsrcs)), byrow=TRUE, ncol=length(q.dsrcs)))

# Decile genomic sizes (in MB)
ys <- setNames(lapply(q.dsrcs, function(ds)
    dcast(x[datasource == ds, .(feature, dataclass, width=width/1e6)], as.integer(feature) ~ dataclass, value.var='width', fun.aggregate=function(x) head(x,1), fill=NA, sep='|')
), q.dsrcs)
pf(ys, col.map=color.map, require.zero=TRUE, lty.map=lty.map)

# _UN_weighted sequencing depth
ys <- setNames(lapply(q.dsrcs, function(ds) {
    sc.dp <- dcast(x[datasource == ds], as.integer(feature) ~ dataclass + group, value.var='unweighted.mean.sc.dp', fun.aggregate=function(x) head(x,1), fill=NA, sep='|')
    #colnames(sc.dp)[-1] <- paste0('sc_', colnames(sc.dp)[-1])
    #bulk.dp <- dcast(x[datasource == ds], as.integer(feature) ~ dataclass + group, value.var='unweighted.mean.bulk.dp', fun.aggregate=function(x) head(x,1), fill=0)
    #colnames(bulk.dp)[-1] <- paste0('bulk_', colnames(bulk.dp)[-1])
    #if (any(sc.dp$feature != bulk.dp$feature)) {
    #   print(sc.dp)
    #   print(bulk.dp)
    #   stop('sc.dp$feature != bulk.dp$feature')
    #}
    #cbind(sc.dp, bulk.dp[,-1])
    # bulk: too much information
    sc.dp
}), q.dsrcs)
pf(ys, col.map=color.map, lty.map=lty.map)


# _UN_weighted allele balance
ys <- setNames(lapply(q.dsrcs, function(ds) dcast(x[datasource == ds], as.integer(feature) ~ dataclass + group, value.var='unweighted.mean.abs.gp.mu', fun.aggregate=function(x) head(x,1), fill=NA, sep='|')), q.dsrcs)
pf(ys, col.map=color.map, lty.map=lty.map)

ys <- setNames(lapply(q.dsrcs, function(ds) dcast(x[datasource == ds], as.integer(feature) ~ dataclass + group, value.var='unweighted.mean.gp.sd', fun.aggregate=function(x) head(x,1), fill=NA, sep='|')), q.dsrcs)
pf(ys, col.map=color.map, lty.map=lty.map)

# Weighted SNV sensitivity
ys <- setNames(lapply(q.dsrcs, function(ds) dcast(x[datasource == ds & muttype == 'snv'], as.integer(feature) ~ dataclass + group, value.var='weighted.sens', sep='|')), q.dsrcs)
pf(ys, col.map=color.map, lty.map=lty.map)

# Weighted indel sensitivity
ys <- setNames(lapply(q.dsrcs, function(ds) dcast(x[datasource == ds & muttype == 'indel'], as.integer(feature) ~ dataclass + group, value.var='weighted.sens', fun.aggregate=function(x) head(x,1), fill=NA, sep='|')), q.dsrcs)
pf(ys, col.map=color.map, lty.map=lty.map, show.xaxis=TRUE)

# Make color legends
# makes a blank plot
blank <- function() { oldmar <- par(mar=c(0,0,0,0)); plot(1, pch=NA, bty='n', yaxt='n', xaxt='n', ylab='', xlab='', main='') }

blank()
this.legend <- color.map[c(1,3:7)]
legend("center", fill=this.legend, legend=names(this.legend), bty='n', border=this.legend)
blank()
blank()
blank()
this.legend <- color.map[c('H3K27ac','H3K36me3','H3K4me1','H3K4me3','H3K9ac')]
legend("center", fill=this.legend, legend=names(this.legend), bty='n', border=this.legend)
blank()
this.legend <- color.map[c('H3K27me3','H3K9me3')]
legend("center", fill=this.legend, legend=names(this.legend), bty='n', border=this.legend)
blank()
this.legend <- c('Brain', 'Non-brain')
legend("center", fill=c('green', 'grey'), legend=this.legend, bty='n', border=c('green', 'grey'), title='GTEx tissue')
blank()
this.legend <- color.map[c('CNS-GBM', 'CNS-Medullo', 'CNS-Oligo',  'CNS-PiloAstro')]
legend("center", fill=this.legend, legend=names(this.legend), bty='n', border=this.legend, title='Tumor type')
dev.off()

if ('snakemake' %in% ls()) {
    sink()
}
