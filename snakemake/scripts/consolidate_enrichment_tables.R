#!/usr/bin/env Rscript

# detect script being run by snakemake
# if so, make a mock commandArgs function
if ('snakemake' %in% ls()) {
    logfile <- snakemake@log[[1]]
    con <- file(logfile, 'w')
    sink(con, type='output')
    sink(con, type='message')

    commandArgs <- function(...) unlist(c(
        snakemake@output['csv'],
        snakemake@params['tag'],
        snakemake@input['snv_indel_pairs']
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 4) {
    cat('each signal file is assumed to come from a single datasource (e.g., scrnaseq).\n')
    cat("`tag` is added to each line to identify the group of cells from which these signals were derived\n")
    cat("This script currently only recognizes signals from datasource in { scrnaseq, scatacseq, encode (dataclass=histone_mark), repliseq }\n")
    stop('usage: consolidate_enrichment_tables.R out.csv tag signal1_snv signal1_indel [ signal2_snv signal2_indel ... signalN_snv signalN_indel ]')
}


out.csv <- args[1]
tag <- args[2]
file.pairs <- args[-(1:2)]
snv.files <- file.pairs[seq(1, length(file.pairs), 2)]
indel.files <- file.pairs[seq(2, length(file.pairs), 2)]

cat(tag, '\n')
cat('snv.files:\n')
print(snv.files)
cat('indel.files:\n')
print(indel.files)

if (file.exists(out.csv))
    stop(paste('output file', out.csv, 'already exists, please delete it first'))

suppressMessages(library(data.table))
suppressMessages(library(mutenrich))


# Read in tables from various genomic covariates. These can have different
# metadata columns, so normalize them to remove metadata columns except:
#   datasource - one of: scrnaseq, scatacseq, repliseq, inactive_histone_mark, active_histone_mark
#   lineclass - for enrichment plots with multiple lines, this column defines each line
# 'tag' adds a column to the table recording what group of cells generated these enrichment signals.
normalize.table <- function(file, tag, ourmuttype=c('snv','indel')) {
    ourmuttype <- match.arg(ourmuttype)

    wrap.fread <- function(...) {
        fread(...)[quantile %in% 1:100][order(as.integer(quantile))][, c('mutsfrom', 'muttype') := list(tag, ourmuttype)][BINSIZE == 1000]
    }

    line1 <- fread(file, nrows=1)
    data.source <- line1[['datasource']]
    if (data.source == 'encode') {
        if (line1[['dataclass']] != 'histone_marks') {
            stop(past('in file', file, 'datasource=', data.source, 'but dataclass=', line1[['dataclass']], '. expected dataclass=histone_marks'))
        }
        # Only one of the 13 brain tissue epigenomes is used
        hst <- wrap.fread(file)[eid == "E073"]
        hst$lineclass <- hst$mark
        hst$datasource <- ifelse(hst$mark %in% c('H3K27me3', 'H3K9me3'), 'inactive_histone_mark', 'active_histone_mark')
        hst$dataclass <- NULL
        hst$signal_type <- NULL
        hst$eid <- NULL
        hst$mark <- NULL
        return(hst)
    } else if (data.source == 'scatacseq') {
        atac <- wrap.fread(file)[libid == 'librarymerged']
        map <- c(OPC='OPCs', astrocyte='Astrocytes', endothelial='Endothelial',
            excitatory_neuron='Excitatory-Neurons', inhibitory_neuron='Inhibitory-Neurons', 
            microglia="Microglia", oligo="Oligodendrocytes")
        atac$celltype <- map[atac$celltype]
        atac$lineclass <- atac$celltype
        atac$libid <- NULL
        atac$sample <- NULL
        atac$celltype <- NULL
        return(atac)
    } else if (data.source == 'scrnaseq') {
        rna <- wrap.fread(file)[donor == 'combined' & selection == 'combined']
        rna$lineclass <- rna$celltype
        rna$dataclass <- NULL
        rna$donor <- NULL
        rna$selection <- NULL
        rna$celltype <- NULL
        return(rna)
    } else if (data.source == 'repliseq') {
        # RepliSeq from ENCODE. 
        repliseq <- wrap.fread(file)
        # have to reverse quantile to make 1=early .. N=late
        max.q <- max(as.integer(repliseq$quantile), na.rm=TRUE)
        min.q <- min(as.integer(repliseq$quantile), na.rm=TRUE)
        repliseq$quantile <- (max.q:min.q)[as.integer(repliseq$quantile)]
        repliseq$lineclass <- paste(repliseq$celltype, repliseq$geoid)
        repliseq$celltype <- NULL
        repliseq$geoid <- NULL
        return(repliseq)
    } else {
        stop(paste('unrecognized datasource ', data.source, ' in file ', file))
    }
}

snv.signals <- rbindlist(lapply(snv.files, normalize.table, tag=tag, ourmuttype='snv'))
indel.signals <- rbindlist(lapply(indel.files, normalize.table, tag=tag, ourmuttype='indel'))

fwrite(rbind(snv.signals, indel.signals), file=out.csv)

if ('snakemake' %in% ls()) {
    sink()
}
