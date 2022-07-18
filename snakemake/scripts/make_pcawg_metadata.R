#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 4)
    stop('usage: make_metadata_helper.R pcawg_sample_sheet.tsv icgc_donor_ids.txt tcga_donor_ids.txt out.txt')

ss <- args[1]
icgc.file <- args[2]
tcga.file <- args[3]
out.csv <- args[4]

if (file.exists(out.csv))
    stop(paste('output file', out.csv, 'already exists, please delete it first'))

library(data.table)

x <- fread(ss)

icgc <- fread(icgc.file)
tcga <- fread(tcga.file)

x[icgc_donor_id %in% icgc[[2]], project := 'icgc']
x[icgc_donor_id %in% tcga[[2]], project := 'tcga']

# Only select tumors with WGS data, not excluded and only retain one record per donor ID
# (some donors have multiple tumor samples, e.g.).
x <- x[library_strategy=='WGS' & grepl('tumour', dcc_specimen_type) & donor_wgs_exclusion_white_gray == 'Whitelist',.SD[1,],by='icgc_donor_id']

y <- rbind(icgc, tcga)
x[y, tumor := Project_Code, on=c('icgc_donor_id'='Donor_ID')]

x <- x[, .(icgc_donor_id, tumor, tumor_code=sapply(strsplit(dcc_project_code, '-'), head, 1), project, maf_processing_error=icgc_donor_id %in% c('DO51542', 'DO35577'), file=paste0('data/pcawg/sample_vcfs/', icgc_donor_id, '.vcf'))]

fwrite(x, file=out.csv)
