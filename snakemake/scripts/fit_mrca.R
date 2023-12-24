#!/usr/bin/env Rscript

# detect script being run by snakemake
# if so, make a mock commandArgs function
if ('snakemake' %in% ls()) {
    logfile <- snakemake@log[[1]]
    con <- file(logfile, 'w')
    sink(con, type='output')
    sink(con, type='message')

    commandArgs <- function(...) unlist(c(
        snakemake@input['object1'],
        snakemake@input['object2'],
        snakemake@input['metadata'],
        snakemake@input['cosmic'],
        snakemake@input['sigb'],
        snakemake@input['aging_model'],
        snakemake@input['aging_model_sigs'],
        snakemake@output['muts'],
        snakemake@output['fits']
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 9) {
    cat("WARNING!! MRCA estimates only apply to oligo-oligo pairs! Only pta_oligo aging rates for total mutation burden and SBS1 are used!\n")
    stop('usage: fit_mrca.R scan2_object1 scan2_object2 sample_metadata.csv cosmic.csv sigb.csv aging_model.csv aging_model_sigs.csv output.muts.csv output.fits.csv')
}


results1.rda <- args[1]
results2.rda <- args[2]
metadata.csv <- args[3]
cosmic.csv <- args[4]
sigb.csv <- args[5]
aging.burden.model.csv <- args[6]
aging.burden.model.sigs.csv <- args[7]
out.muts.csv <- args[8]
out.fits.csv <- args[9]

for (f in c(out.muts.csv, out.fits.csv)) {
    if (file.exists(f))
        stop(paste('output file', f, 'already exists, please delete it first'))
}

suppressMessages(library(scan2))
suppressMessages(library(pracma))

combine.dfs <- function(s1, s2) {
    # just rename some columns to be unambiguous after cbind()
    s1$pass1 <- s1$pass
    s1$af1 <- s1$af
    s1$dp1 <- s1$dp
    s1$alt1 <- ifelse(s1$dp == 0, 0, s1$af*s1$dp)

    s2$pass2 <- s2$pass
    s2$af2 <- s2$af
    s2$dp2 <- s2$dp
    s2$alt2 <- ifelse(s2$dp == 0, 0, s2$af*s2$dp)

    d1 <- cbind(s1, s2)
    d1$shared <- (d1$pass1 | d1$pass2) & d1$alt1 > 1 & d1$alt2 > 1
    d1$private <- !d1$shared & (d1$pass1 | d1$pass2) &
        (d1$alt1 == 0 | d1$alt2 == 0) &   # added may 25 2022
        (d1$dp1 > 5 & d1$dp2 > 5)         # added may 25 2022
    # anything not shared or private, but must be passed somewhere
    d1$indeterminate <- (d1$pass1 | d1$pass2) & !d1$shared & !d1$private  
    # labels are just convenient for tabling and grouping
    d1$classification <- ifelse(d1$shared, 'shared', ifelse(d1$private, 'private', ifelse(d1$indeterminate, 'indeterminate', 'none')))
    d1
}

aging.burden.model <- fread(aging.burden.model.csv)[Group == 'pta_oligo']
aging.burden.model.sigs <- fread(aging.burden.model.sigs.csv)[Group == 'pta_oligo' & Sig =='SBS1']

metadata <- fread(metadata.csv)

cosmic <- fread(cosmic.csv)
pta.artifact <- get(data(snv.artifact.signature.v3))   # loads `snv.artifact.signature.v3`
sigb <- fread(sigb.csv)$B

cat('loading SCAN2 object 1..\n')
r1 <- get(load(results1.rda))
amptype1 <- metadata[sample == r1@single.cell]$amp
cat('loading SCAN2 object 2..\n')
r2 <- get(load(results2.rda))
amptype2 <- metadata[sample == r2@single.cell]$amp
cat('memory usage:\n')
gc()

x <- combine.dfs(r1@gatk[somatic.candidate == TRUE & muttype == 'snv'],
                 r2@gatk[somatic.candidate == TRUE & muttype == 'snv'])
#muttab <- x[shared == TRUE | private == TRUE]
muttab <- x[classification != 'none']

fwrite(muttab, file=out.muts.csv)


# Use leave-one-out hSNPs to quantify misclassification of shared as private.
# All hSNPs are should be shared.
hsnptab <- combine.dfs(r1@gatk[resampled.training.site == TRUE & muttype == 'snv'][, pass := training.pass],
                       r2@gatk[resampled.training.site == TRUE & muttype == 'snv'][, pass := training.pass])

# IMPORTANT: because `cosmic` contains SigB for MDA samples, there is no
# need to further correct for sig B when using get.sbs1().
get.sbs1 <- function(tab) {
    setNames(lsqnonneg(as.matrix(cosmic[,-1]), as.vector(table(sbs96(tab$mutsig))))$x,
        colnames(cosmic)[-1])['SBS1']
}

# WARNING!!! muttab and hsnptab MUST be subsampled to pass==TRUE for the
# cell corresponding to `object`!!
# override n.shared and n.private to do SBS1 modeling
# this.cosmic is the global cosmic plus either the PTA artifact signature or MDA artifact signature (sigB)
f <- function(object, amptype, this.cosmic, muttab, hsnptab, model.intercept, model.intercept.lb, model.intercept.ub, model.rate, model.rate.lb, model.rate.ub, remove.sigb=FALSE,
    n.shared.fn=nrow, n.private.fn=nrow)
{
    sample <- object@single.cell
    burden <- mutburden(object)['snv']
    n.shared <- n.shared.fn(muttab[shared == TRUE])
    n.private <- n.private.fn(muttab[private == TRUE])

    # Quantity is: number of shareds misclassified as a function of total shared count.
    # This doesn't matter much for PTA, but does for MDA.
    # When doing SBS1, no need to estimate SBS1 level and/or restrict hSNPs to SBS1-like muts
    print(table(private=hsnptab$private, shared=hsnptab$shared))

    # this doesn't make sense, but also we need to account for indeterminate
    #shared.as.private.rate <- sum(hsnptab$private == TRUE, na.rm=TRUE) / sum(hsnptab$shared == TRUE, na.rm=TRUE)
    # all germline hSNPs are shared mutations. anything not classified as shared is considered
    # a misclassification (whether it is indeterminate or private).
    shared.misclassification.rate <- hsnptab[, mean(classification != 'shared')]

    # Shift the estimated number of shared mutations accidentally misclassified as private
    # from the private count to the shared count.
    # 
    # T: unobserved.  True number of shared mutations between cells (unobserved)
    # S: observed.    Number of shared mutations detected by the heuristic (observed)
    # f: inferred.    Fraction of truly shared mutations T misclassified by heuristic as private or indeterminate
    # fT:             Number of misclassified shared SNVs as private or indet.
    # 
    # Total number of shared mutations is the sum of the ones we detect as shared and
    # those misclassified as private.
    #
    #   T = S + fT
    #   T - fT = S
    #   T(1 - f) = S
    #   T = S/(1 - f)
    n.shared.adj <- n.shared / (1 - shared.misclassification.rate)
    n.private.adj <- n.private - (n.shared.adj - n.shared)

    # For MDA samples, private mutations will have signature B artifacts.
    # shared (mostly) are not affected by sig B.
    sigb <- 0
    if (remove.sigb) {
        fit <- setNames(lsqnonneg(as.matrix(this.cosmic[,-1]), as.vector(table(sbs96(muttab[private == TRUE]$mutsig))))$x, colnames(this.cosmic)[-1])
        sigb <- fit['SigB'] / sum(fit)
    }
    n.private.adj2 <- n.private.adj * (1 - sigb)
    burden.adj <- burden * (1 - sigb)
    mut.to.burden.rate <- burden.adj / nrow(object@gatk[muttype == 'snv' & pass == TRUE])

    # Extrapolate rates to whole genome burden. This matches the rates used
    # in the aging models.
    n.shared.final <- n.shared.adj * mut.to.burden.rate
    n.private.final <- n.private.adj2 * mut.to.burden.rate

    # New method: use the variability in the aging model's intercept/yearly rate to
    # create an interval on the number of shared years. This method is probably better
    # than the original count-forward, count-backward method because shared mutations
    # are (1) far more likely to not be artifactual and (2) it is very easy to model
    # sensitivity since all hSNPs are also shared and can be used to accurately
    # determine what fraction of shared sites should be recovered.
    #
    # All estimates are of the form:
    #
    #    years  =  (# of sensitivity adjusted shared mutations - INTERCEPT) / RATE
    #
    # where INTERCEPT and RATE each take on the point estimate from the linear model
    # and the 95% CI upper and lower bounds.

    # Get all possible combinations of upper and lower bounds. The min and max values
    # don't correspond to min and max values for the final range because the values
    # (years) can be negative. year=0 means birth, so there are still 9 months
    # preceding that.
    grid <- expand.grid(c(model.intercept.lb, model.intercept.ub),
                c(model.rate.lb, model.rate.ub))
    ests <- apply(grid, 1, function(row) (n.shared.final - row[1]) / row[2])
    cat("all estimates using 95% CI bounds:\n")
    print(ests)

    # As explained above, these are not actually the min/max values because n.shared.final-intercept can be <0
    # Maximum (=upper bounds) INTERCEPT and RATE => minimum value for years (=lower bound)
    #years.shared.lb <- (n.shared.final - model.intercept.ub) / model.rate.ub
    # Minimum (=lower bounds) INTERCEPT and RATE => maximum value for years (=upper bound)
    #ears.shared.ub <- (n.shared.final - model.intercept.lb) / model.rate.lb
    years.shared.lb <- min(ests)
    years.shared.ub <- max(ests)
    years.shared <- (n.shared.final - model.intercept) / model.rate

    # no longer relevant, but kept around for posterity
    years.private <- (n.private.final) / model.rate
    years.total <- (n.shared.final + n.private.final - model.intercept) / model.rate

    data.frame(sample=sample, amptype=amptype, sigb=sigb, shared.misclassification.rate=shared.misclassification.rate,
        burden=burden, burden.adj=burden.adj, mut.to.burden.rate=mut.to.burden.rate,
        n.shared=n.shared, n.shared.adj=n.shared.adj, n.shared.final=n.shared.final,
        n.private=n.private, n.private.adj=n.private.adj, n.private.adj2=n.private.adj2, n.private.final=n.private.final,
        model.intercept=model.intercept,
        model.intercept.lb=model.intercept.lb,
        model.intercept.ub=model.intercept.ub,
        model.rate=model.rate,
        model.rate.lb=model.rate.lb,
        model.rate.ub=model.rate.ub,
        years.shared=years.shared,
        years.shared.lb=years.shared.lb,
        years.shared.ub=years.shared.ub,
        years.private=years.private,   # these are no longer relevant
        years.total=years.total)       #
}

# Date MRCA based on total mutation burden
# all.shared is the union of shared sites. the n.shared values from f()
# are the number of shared AND called sites in a specific sample. this is
# the relevant value for extrapolating since sensitivity is relatively well
# approximated for single-sample somatic calling with SCAN2.

pta.cosmic <- cbind(cosmic, PTA=pta.artifact)
mda.cosmic <- cbind(cosmic, SigB=sigb)
cosmics <- list(PTA=pta.cosmic, MDA=mda.cosmic)

total.mut <- data.frame(model='Total burden', all.shared=sum(muttab$shared),
    rbind(f(object=r1,
            amptype=amptype1,
            this.cosmic=cosmics[[amptype1]],
            muttab=muttab[pass1 == TRUE],
            hsnptab=hsnptab[pass1 == TRUE],
            model.intercept=aging.burden.model[Variable=='(Intercept)', Estimate],
            model.intercept.lb=aging.burden.model[Variable=='(Intercept)', `95% CI lower`],
            model.intercept.ub=aging.burden.model[Variable=='(Intercept)', `95% CI upper`],
            model.rate=aging.burden.model[Variable=='age', Estimate],
            model.rate.lb=aging.burden.model[Variable=='age', `95% CI lower`],
            model.rate.ub=aging.burden.model[Variable=='age', `95% CI upper`],
            remove.sigb=amptype1 == 'MDA'),
        f(object=r2,
            amptype=amptype2,
            this.cosmic=cosmics[[amptype2]],
            muttab=muttab[pass2 == TRUE],
            hsnptab=hsnptab[pass2 == TRUE],
            model.intercept=aging.burden.model[Variable=='(Intercept)', Estimate],
            model.intercept.lb=aging.burden.model[Variable=='(Intercept)', `95% CI lower`],
            model.intercept.ub=aging.burden.model[Variable=='(Intercept)', `95% CI upper`],
            model.rate=aging.burden.model[Variable=='age', Estimate],
            model.rate.lb=aging.burden.model[Variable=='age', `95% CI lower`],
            model.rate.ub=aging.burden.model[Variable=='age', `95% CI upper`],
            remove.sigb=amptype2 == 'MDA')))

# Date MRCA based on SBS1
sbs1.mut <- data.frame(model='SBS1 burden', all.shared=NA,
    rbind(f(object=r1,
            amptype=amptype1,
            this.cosmic=cosmics[[amptype1]],
            muttab=muttab[pass1 == TRUE],
            hsnptab=hsnptab[pass1 == TRUE],
            model.intercept=aging.burden.model.sigs[Variable == '(Intercept)', Estimate],
            model.intercept.lb=aging.burden.model.sigs[Variable=='(Intercept)', `95% CI lower`],
            model.intercept.ub=aging.burden.model.sigs[Variable=='(Intercept)', `95% CI upper`],
            model.rate=aging.burden.model.sigs[Variable == 'age', Estimate],
            model.rate.lb=aging.burden.model.sigs[Variable=='age', `95% CI lower`],
            model.rate.ub=aging.burden.model.sigs[Variable=='age', `95% CI upper`],
            n.shared.fn=get.sbs1, n.private.fn=get.sbs1),
        f(object=r2,
            amptype=amptype2,
            this.cosmic=cosmics[[amptype2]],
            muttab=muttab[pass2 == TRUE],
            hsnptab=hsnptab[pass2 == TRUE],
            model.intercept=aging.burden.model.sigs[Variable == '(Intercept)', Estimate],
            model.intercept.lb=aging.burden.model.sigs[Variable=='(Intercept)', `95% CI lower`],
            model.intercept.ub=aging.burden.model.sigs[Variable=='(Intercept)', `95% CI upper`],
            model.rate=aging.burden.model.sigs[Variable == 'age', Estimate],
            model.rate.lb=aging.burden.model.sigs[Variable=='age', `95% CI lower`],
            model.rate.ub=aging.burden.model.sigs[Variable=='age', `95% CI upper`],
            n.shared.fn=get.sbs1, n.private.fn=get.sbs1)))

results <- rbind(total.mut, sbs1.mut)
rownames(results) <- NULL

fwrite(results, file=out.fits.csv)

if ('snakemake' %in% ls()) {
    sink()
}
