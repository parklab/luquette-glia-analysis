#!/usr/bin/env Rscript

# detect script being run by snakemake
# if so, make a mock commandArgs function
if ('snakemake' %in% ls()) {
    logfile <- snakemake@log[[1]]
    con <- file(logfile, 'w')
    sink(con, type='output')
    sink(con, type='message')

    commandArgs <- function(...) unlist(c(
        snakemake@input['rda'],
        snakemake@input['sens'],
        snakemake@params['group'],
        snakemake@params['muttype'],
        snakemake@params['passtype'],
        snakemake@output['full_rda'],
        snakemake@output['summary_rda']
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}


args <- commandArgs(trailingOnly=T)

if (length(args) != 7)
    stop("usage: apply_sensitivity_correction.R in.full.enrich.rda in.sens.table group {snv|indel} {A|AB} out.full.rda out.summary.rda")

e.full.rda <- args[1]
sens.file <- args[2]
this.group <- args[3]
this.muttype <- args[4]
this.passtype <- args[5]
out.full.rda <- args[6]
out.summary.rda <- args[7]

for (f in c(out.full.rda, out.summary.rda)) {
    if (file.exists(f))
        stop(sprintf("output file %s already exists, please delete it first", f))
}

if (this.passtype == 'A') {
    this.passtype <- 'pass'
} else if (this.passtype == 'AB') {
    this.passtype <- 'rescue'
} else {
    stop(paste('unrecognized passtype', this.passtype, 'must be A or AB (case sensitive)'))
}


suppressMessages(library(scan2))
suppressMessages(library(mutenrich))

apply.correction <- function(e, emeta, sens) {
    require(mutenrich)

    # allow for joining and so on, but return a modified version of the original
    # data.frame at the end of this function to maintain complete compatibility
    # with older code.
    emeta.dt <- as.data.table(emeta)

    # enrichment analyses are always run on signals from the same datasource.
    dsrc <- emeta.dt$datasource[1]
    dclass <- emeta.dt$dataclass[1]
    # repliseq has no dataclass
    if (length(dclass) == 0)
        dclass <- ''

    # Signature enrichment is detectable by the `signature.catalog` stored in
    # the enrichment object.
    sigenrich <- ('signature.catalog' %in% names(e$edata))
    if (sigenrich) {
        sigcat <- e$edata$signature.catalog
        cat("detected signature enrichment using catalog of", ncol(sigcat), "signatures:\n")
        print(colnames(sigcat))
    }

    # Determine the subset of enrichment analyses that we want to correct with our sensitivity table.
    # The reason this is not all enrichment analyses in the object is that sensitivity metrics were
    # only collected for the final, reported signals, which are a subset of those tested.
    this.datasource <- dsrc
    this.dataclass <- dclass
    emeta.dt[, signal := e$feature.set]
    if (dsrc == 'scatacseq') {
        this.subset <- emeta.dt[, .(subset=libid == 'librarymerged' & sample == 'merged')]$subset
        this.dataclass <- 'celltype'
    } else if (dsrc == 'scrnaseq') {
        this.subset <- emeta.dt[, .(subset=selection == 'combined' & donor == 'combined')]$subset
        this.dataclass <- 'celltype'
    } else if (dsrc == 'encode' & dclass == 'histone_marks') {
        this.subset <- emeta.dt[, .(subset=eid == 'E073')]$subset
        this.dataclass <- 'mark'
        this.datasource <- c('active_histone_mark', 'inactive_histone_mark')
        emeta.dt[, datasource := ifelse(mark %in% c('H3K9me3', 'H3K27me3'), 'inactive_histone_mark', 'active_histone_mark')]
    } else if (dsrc == 'repliseq') {
        this.subset <- emeta.dt[, .(subset=datasource == 'repliseq')]$subset  # this just gets everything
        this.dataclass <- 'celltype'
    } else if (dsrc == 'gtex') {
        this.subset <- emeta.dt[, .(subset=datasource == 'gtex')]$subset  # this just gets everything
        # datasource=gtex can either be gtex transcription or gtex tx/utx regions.
        # these two types are never mixed, so just check the first element to see
        # which this is.
        # N.B., gtex expression doesn't have a dataclass column, so can't check that
        if ('dataclass' %in% colnames(emeta.dt[this.subset,])) {
            if (emeta.dt[this.subset,]$dataclass[1] == 'tx_or_utx') {
                this.dataclass <- 'dataclass'
            }
        } else {
            this.dataclass <- 'tissue'
        }
    } else if (dsrc == 'gencode') {
        this.subset <- emeta.dt[, .(subset=datasource == 'gencode')]$subset  # this just gets everything
        this.dataclass <- 'dataclass'
    } else if (dsrc == 'gencode_simplified') {
        this.subset <- emeta.dt[, .(subset=datasource == 'gencode_simplified')]$subset  # this just gets everything
        this.dataclass <- 'dataclass'
    } else if (dsrc == 'nott_enhprom') {
        this.subset <- emeta.dt[, .(subset=datasource == 'nott_enhprom')]$subset  # this just gets everything
        this.dataclass <- 'dataclass'
    } else if (dsrc == 'roadmap') {
        this.subset <- emeta.dt[, .(subset=datasource == 'roadmap' & eid == 'E073')]$subset
        this.dataclass <- 'dataclass'
    } else {
        stop(paste('unknown datasource', dsrc, 'and dataclass', dclass))
    }

    sens <- sens[datasource %in% this.datasource]
    if (nrow(sens) == 0)
        stop(paste('datasource', this.datasource, 'in enrichment object is not found in sensitivity table'))


    emeta <- emeta[this.subset,]  # also subset the original data.frame for the final return value
    emeta.dt <- emeta.dt[this.subset]
    emeta.dt$dataclass <- emeta.dt[[this.dataclass]]
    # mutenrich::subset() below does not update this. it was intended to subset by features
    # within one signal, not subsetting at the signal level.
    attr(e$edata$gbed, 'feature.set') <- attr(e$edata$gbed, 'feature.set')[this.subset]

    # Now map each (signal,feature) pair to a corrected observation count by
    # adjusting for the estimated sensitivity of the (signal,feature) region.
    #
    # ASSUMPTION: QBEDs and regular BEDs are never mixed. This is true for now: different
    # scripts handle each (qbedenrich.R, bedenrich.R).
    # If emeta.dt contains QUANTILES, MIN_COVERAGE and BINSIZE, assume it is a QBED.
    if (all(c('QUANTILES', 'MIN_COVERAGE', 'BINSIZE') %in% colnames(emeta.dt))) {
        cat('detected QBED: joining on n.quantiles and binsize\n')
        sens[, n.quantiles := as.character(n.quantiles)]
        sens[, binsize := as.character(binsize)]
        sens <- sens[emeta.dt,,on=.(datasource,dataclass,n.quantiles=QUANTILES,binsize=BINSIZE)]
    } else {
        cat('detected BED\n')
        sens <- sens[emeta.dt,,on=.(datasource,dataclass)]
    }

    if (sigenrich) {
        # The most convenient way to handle this is to duplicate the row in the
        # table for each feature (which usually means a quantile of a QBED signal)
        # for each signature in the catalog. I.e., the same sensitivity estimate
        # is used for all signatures.
        #
        # This is valid because signature enrichment is never performed on SCAN2
        # rescue mutations, only VAF-based mutations, which should not be biased
        # by mutation signatures.
        #
        # In addition to the signatures in the catalog, there are some measures of  
        # goodness of fit of the signature reconstruction.
        all.feats <- c(colnames(sigcat), 'resid.norm', 'norm', 'pct.resid')
        sens <- rbindlist(lapply(all.feats, function(signature.feat) {
            x <- copy(sens)  # data.tables are weird about copy mechanics
            x[, feature := paste0(feature, '|||', signature.feat)]
            # also need to store the signature separately for easy grouping
            x[, signature.feat := signature.feat]
            x   
        }))

        # The difference between signature enrichment and regular enrichment is
        # that sensitivity correction is applied within each signature (which is a
        # subset of the total mutation calls) and that we need to scale the value of
        # `obs`, not n.calls.  `obs` is signature exposure, not a count of mutations.
        # however, it is still correct to use n.calls for calculating the scaling
        # factor because n.calls is how many mutations are permuted around the
        # genome in the permutation analysis.
        #
        # N.B. the sensitivity table only knows about the counts of mutations used
        # for permutation analysis. It does not know how much of each signature was
        # fit to each bin. So, unlike the non-sigenrich analysis, we do not have the
        # necessary info to compute the corrected final value in correction.map.
        # Instead, store the scaling factor and apply it later.
        correction.map <-
            sens[, .(sigfeat=paste0(signal, '|||', feature),
                    # remember: don't confuse this value with the signature exposure level in `obs`
                    n.calls=n.calls,   
                    weighted.sens=ifelse(is.na(weighted.sens), 1, weighted.sens),
                    weighted.sum=sum(n.calls/ifelse(is.na(weighted.sens), 1, weighted.sens)),
                    # this is the SCALING FACTOR, not the final corrected value, so n.calls is replaced by '1'
                    corr.n.calls=1/ifelse(is.na(weighted.sens), 1, weighted.sens) * sum(n.calls)/sum(n.calls/ifelse(is.na(weighted.sens), 1, weighted.sens))),
                by=.(dataclass, signature.feat)]

            # DO NOT scale the norms of the matrix or residual
            # Of course, these values are no longer correct, but scaling them makes little sense.
            # It might be better to just set them to NA to prevent mis-interpretation.
            sens[signature.feat %in% c('norm', 'pct.resid', 'resid.norm'), corr.n.calls := 1]
    } else {
        correction.map <-
            sens[, .(sigfeat=paste0(signal, '|||', feature),
                    # record these values for sanity checks and debugging
                    n.calls=n.calls,
                    weighted.sens=ifelse(is.na(weighted.sens), 1, weighted.sens),
                    weighted.sum=sum(n.calls/ifelse(is.na(weighted.sens), 1, weighted.sens)),
                    corr.n.calls=n.calls / ifelse(is.na(weighted.sens), 1, weighted.sens) * sum(n.calls)/sum(n.calls/ifelse(is.na(weighted.sens), 1, weighted.sens))),
                by=dataclass]
    }
    setkey(correction.map, sigfeat)
    
    # Select the (signal,feature) pairs present in the sensitivity table from the enrichment object.
    # The sensitivity table should be a subset of the enrichment object--some extra enrichment signals
    # were computed that were never used later.
    e <- mutenrich::subset(e, fts=correction.map$sigfeat)

    if (sigenrich) {
        e$real.obs.uncorrected <- e$real.obs  # can't use correction.map like the non-sigenrich case
        # "corr.n.calls" is the scaling factor for sigenrich analysis
        e$real.obs <- setNames(e$real.obs * correction.map[names(e$real.obs)]$corr.n.calls, names(e$real.obs))

        # bootstrapping was not applied to signature enrichment because label-shuffling doesn't work there.
        cat('sigenrich: skipping bootstrap\n')
    } else {
        # Sanity check: ensure the observed counts in the enrichment object (e$real.obs)
        # are identical to the call counts in the sensitivity table (correction.map$n.calls).
        # This does not apply to signature enrichment, in which "counts" are the least squares
        # weights of the various signatures.
        if (!all(correction.map[names(e$real.obs)]$n.calls == e$real.obs)) {
            stop('n.calls in sensitivity table do not all match real.obs in enrichment object')
        }
        e$real.obs.uncorrected <- setNames(correction.map[names(e$real.obs)]$n.calls, names(e$real.obs)) # store for posterity
        # round(): the observation values must be integers for downstream bootstrapping
        e$real.obs <- setNames(round(correction.map[names(e$real.obs)]$corr.n.calls), names(e$real.obs))

        # Apply bootstrapping
        n.boot <- nrow(e$boots)  # all signals/features receive same number of bootstraps
        new.boots <- corrected.bootstrap(e=e, edata=e$edata, n.boot=n.boot)
        # not all features are in this output because of the subset() above. the sensitivity
        # table does not have, e.g., the "outside" feature for QBEDs since this exclusively
        # refers to chrX and chrY. ("outside" is in the sensitivity table for some BEDs)
        boots.cols.to.keep <- names(e$real.obs)
        e$boots <- new.boots[,..boots.cols.to.keep]
    }

print(correction.map)
print(e$real.obs)
    es <- setNames(lapply(feature.set(edata=e$edata), function(feat.name) {
        this.e <- get.feat(e, feat.name)
        esummary(this.e)
    }), feature.set(edata=e$edata))

    list(e=e, es=es, emeta=emeta, sens=sens, correction.map=correction.map, sigenrich=sigenrich)
}

# For bootstrapping, need to use the sensitivity corrected obs values, not the original
# values. Since the corrected values only represent a sensitivity-corrected *count* of
# mutations, there is no actual mutation set corresponding to the corrected values.
# So this map function just returns a dummy list of labels with the appropriate frequency.
#
# This is an exact copy of the mutenrich::bootstrap function except where noted (real.labels).
#
# Original bootstrap() is function(grl,...) rather than function(e,...). We need 'e'
# for the corrected obs values.
corrected.bootstrap <- function(e, edata, n.boot, verbose=1) {
    if (verbose > 0) cat(sprintf('Bootstrapping %d samples: [', n.boot))
    if (verbose == 2) { cat('gc3a\n'); print(gc()) }

    ret <- do.call(cbind, lapply(feature.set(edata=edata), function(feat.name) {
        # real.labels used to be taken from the observed mutations. now that we only
        # know the corrected counts, the label set here should just be a vector of
        # features with the corrected frequencies.
        #real.labels <- edata$map.fn(GenomicRanges::GRangesList(gmuts),
            #edata=edata, feat.name=feat.name)$to
        this.e <- get.feat(en=e, feat.name=feat.name, strip.feat.name=TRUE)
        real.labels <- factor(
            rep(names(this.e$real.obs), this.e$real.obs),
            # important: use the same factor encoding
            levels=levels(features(edata=edata, feat.name=feat.name))
        )

        # the combined dataset has to have an ID called 'perm.id' for
        # map.and.count to group the mapping prior to counting.
        if (verbose == 2) cat('sampling', length(real.labels)*n.boot, 'for bootstrapping\n')
        boot.obs <- data.table::data.table(t(sapply(1:n.boot, function(i) {
            if (i %% (n.boot/10) == 0 & verbose==1) cat('.')
            s <- sample(real.labels, size=length(real.labels), replace=TRUE)
            edata$count.fn(s, edata, feat.name)
        })))
        if (verbose == 2) { cat('gc3b\n'); print(gc()) }
        boot.obs
    }))
    if (verbose > 0) cat('] 100%')
    ret
}

sens.tab <- fread(sens.file)

load(e.full.rda)


# Save for comparison
e.uncorrected <- e
es.uncorrected <- es
emeta.uncorrected <- emeta

corrected.output <- apply.correction(e=e, emeta=emeta,
    sens=sens.tab[group == this.group & muttype == this.muttype & passtype == this.passtype])
e <- corrected.output$e
es <- corrected.output$es
emeta <- corrected.output$emeta
correction.map <- corrected.output$correction.map
sigenrich <- corrected.output$sigenrich
sens <- corrected.output$sens

save(e, es, emeta, e.uncorrected, es.uncorrected, emeta.uncorrected, sens, correction.map, sigenrich, file=out.full.rda)
save(es, emeta, es.uncorrected, emeta.uncorrected, sens, correction.map, sigenrich, file=out.summary.rda)

if ('snakemake' %in% ls()) {
    sink()
}
