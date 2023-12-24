This directory is only for external data files that cannot be automatically downloaded.


NanoSeq and META-CS data might be automatically downloadable--have not checked.

NeuronDuplexSequencingComparison.xlsx was manually constructed from the Abascal NanoSeq
and Xing META-CS data.


From_Javier:
    Excel spreadsheet from Javier. Contents of each tab were manually pasted into text files.
    Panel 3 - final 20 lines are GENENAME_MUTTYPE where MUTTYPE is something like "CAG" for
        (I guess) CAG expansion or "-Gain" for a copy gain (again, I guess). The _MUTTYPE
        suffixes were manually deleted (to retain the gene).
    Finally, each file was sorted and uniqed to remove duplicate genes.

    These 4 files were combined into a single BED file with record names reflecting the panel.
    Gene names were matched to our gene model by joining via data.table

        library(data.table)
        g <- fread("../data/gtex/gencode.v26lift37.annotation.GTEX_COLLAPSED.genes_only.bed")
        ps <- lapply(list.files(pattern='.*.uniq_sort.txt'), function(f) { x <- fread(f,header=F); colnames(x) <- 'gene'; x })      
        ps2 <- lapply(ps, function(p) g[p,,on='gene'])     
        ps3 <- lapply(ps2, function(p) p[!is.na(chr)])
        ps4 <- lapply(1:4, function(i) ps3[[i]][, panel := paste0('panel', i)])
        fwrite(rbindlist(ps4), file='combined_4_panels.bed', sep='\t', col.names=FALSE)

    The resulting matches were good, especially given that we do not include chrX or chrY in
    our gene model:

        panel        1    2    3    4  (in order of list.files())
        not found    6   30   33    7
        Total      118  311  414   98
