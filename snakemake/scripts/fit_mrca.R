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
        snakemake@input['cosmic'],
        snakemake@input['sigb'],
        snakemake@input['aging_model'],
        snakemake@input['aging_model_sigs'],
        snakemake@params['amptype'],
        snakemake@output['muts'],
        snakemake@output['fits']
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 9) {
    stop('usage: fit_mrca.R scan2_object1 scan2_object2 cosmic.csv sigb.csv aging_model.csv aging_model_sigs.csv {mda|pta} output.muts.csv output.fits.csv')
}


results1.rda <- args[1]
results2.rda <- args[2]
cosmic.csv <- args[3]
sigb.csv <- args[4]
aging.burden.model.csv <- args[5]
aging.burden.model.sigs.csv <- args[6]
amptype <- args[7]
out.muts.csv <- args[8]
out.fits.csv <- args[9]

if (amptype != 'mda' & amptype != 'pta')
    stop("amptype must be either 'mda' or 'pta', case sensitive")

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
    d1
}

aging.burden.model <- fread(aging.burden.model.csv)[CellType == 'Oligo']
aging.burden.model.sigs <- fread(aging.burden.model.sigs.csv)[CellType == 'Oligo']

cosmic <- fread(cosmic.csv)
if (amptype == 'pta') {
    data(snv.artifact.signature.v3)  
    cosmic$PTA <- snv.artifact.signature.v3
} else if (amptype == 'mda') {
    lod <- fread(sigb.csv)
    cosmic$SigB <- lod$B
} else {
    stop("amptype must be either 'mda' or 'pta', case sensitive")
}

cat('loading SCAN2 object 1..\n')
r1 <- get(load(results1.rda))
cat('loading SCAN2 object 2..\n')
r2 <- get(load(results2.rda))
cat('memory usage:\n')
gc()

x <- combine.dfs(r1@gatk[somatic.candidate == TRUE & muttype == 'snv'],
                 r2@gatk[somatic.candidate == TRUE & muttype == 'snv'])
muttab <- x[shared == TRUE | private == TRUE]

fwrite(muttab, file=out.muts.csv)


# Use leave-one-out hSNPs to quantify misclassification of shared as private.
# All hSNPs are should be shared.
hsnptab <- combine.dfs(r1@gatk[resampled.training.site == TRUE & muttype == 'snv'][, pass := resampled.training.pass],
                       r2@gatk[resampled.training.site == TRUE & muttype == 'snv'][, pass := resampled.training.pass])

# IMPORTANT: because `cosmic` contains SigB for MDA samples, there is no
# need to further correct for sig B when using get.sbs1().
get.sbs1 <- function(tab) {
    setNames(lsqnonneg(as.matrix(cosmic[,-1]), as.vector(table(sbs96(tab$mutsig))))$x,
        colnames(cosmic)[-1])['SBS1']
}

# WARNING!!! muttab and hsnptab MUST be subsampled to pass==TRUE for the
# cell corresponding to `object`!!
# override n.shared and n.private to do SBS1 modeling
f <- function(object, muttab, hsnptab, model.intercept, model.rate, remove.sigb=FALSE,
    n.shared.fn=nrow, n.private.fn=nrow)
{
    sample <- object@single.cell
    burden <- object@mutburden$snv[2,]$burden
    n.shared <- n.shared.fn(muttab[shared == TRUE])
    n.private <- n.private.fn(muttab[private == TRUE])

    # Quantity is: number of shareds misclassified as a function of total shared count.
    # This doesn't matter much for PTA, but does for MDA.
    # When doing SBS1, no need to estimate SBS1 level and/or restrict hSNPs to SBS1-like muts
    print(table(private=hsnptab$private, shared=hsnptab$shared))
    shared.as.private.rate <- sum(hsnptab$private == TRUE, na.rm=TRUE) / sum(hsnptab$shared == TRUE, na.rm=TRUE)

    # Shift the estimated number of shared mutations accidentally misclassified as private
    # from the private count to the shared count.
    n.shared.adj <- n.shared * (1 + shared.as.private.rate)
    n.private.adj <- n.private - (n.shared.adj - n.shared)

    # For MDA samples, private mutations will have signature B artifacts.
    # shared (mostly) are not affected by sig B.
    sigb <- 0
    if (remove.sigb) {
        fit <- setNames(lsqnonneg(as.matrix(cosmic[,-1]), as.vector(table(sbs96(muttab[private == TRUE]$mutsig))))$x, colnames(cosmic)[-1])
        sigb <- fit['SigB'] / sum(fit)
    }
    n.private.adj2 <- n.private.adj * (1 - sigb)
    burden.adj <- burden * (1 - sigb)
    mut.to.burden.rate <- burden.adj / nrow(object@gatk[muttype == 'snv' & pass == TRUE])

    # Extrapolate rates to whole genome burden. This matches the rates used
    # in the aging models.
    n.shared.final <- n.shared * mut.to.burden.rate
    n.private.final <- n.private * mut.to.burden.rate

    years.shared <- (n.shared.final - model.intercept) / model.rate
    years.private <- (n.private.final) / model.rate
    years.total <- (n.shared.final + n.private.final - model.intercept) / model.rate

    data.frame(sample=sample, amptype=amptype, sigb=sigb, shared.as.private.rate=shared.as.private.rate,
        burden=burden, burden.adj=burden.adj, mut.to.burden.rate=mut.to.burden.rate,
        n.shared=n.shared, n.shared.adj=n.shared.adj, n.shared.final=n.shared.final,
        n.private=n.private, n.private.adj=n.private.adj, n.private.adj2=n.private.adj2, n.private.final=n.private.final,
        model.intercept=model.intercept, model.rate=model.rate,
        years.shared=years.shared, years.private=years.private, years.total=years.total)
}

# Date MRCA based on total mutation burden
# all.shared is the union of shared sites. the n.shared values from f()
# are the number of shared AND called sites in a specific sample. this is
# the relevant value for extrapolating since sensitivity is relatively well
# approximated for single-sample somatic calling with SCAN2.
total.mut <- data.frame(model='Total burden', all.shared=sum(muttab$shared),
    rbind(f(r1, muttab[pass1 == TRUE], hsnptab[pass1 == TRUE],
            model.intercept=aging.burden.model[Variable=='Intercept', Estimate],
            model.rate=aging.burden.model[Variable=='Age', Estimate],
            remove.sigb=amptype == 'mda'),
        f(r2, muttab[pass2 == TRUE], hsnptab[pass2 == TRUE],
            model.intercept=aging.burden.model[Variable=='Intercept', Estimate],
            model.rate=aging.burden.model[Variable=='Age', Estimate],
            remove.sigb=amptype == 'mda')))

# Date MRCA based on SBS1
sbs1.mut <- data.frame(model='SBS1 burden', all.shared=NA,
    rbind(f(r1, muttab[pass1 == TRUE], hsnptab[pass1 == TRUE],
            model.intercept=aging.burden.model.sigs[Sig=='SBS1', AgeIntercept],
            model.rate=aging.burden.model.sigs[Sig=='SBS1',AgeRate],
            n.shared.fn=get.sbs1, n.private.fn=get.sbs1),
        f(r2, muttab[pass2 == TRUE], hsnptab[pass2 == TRUE],
            model.intercept=aging.burden.model.sigs[Sig=='SBS1', AgeIntercept],
            model.rate=aging.burden.model.sigs[Sig=='SBS1',AgeRate],
            n.shared.fn=get.sbs1, n.private.fn=get.sbs1)))

results <- rbind(total.mut, sbs1.mut)
rownames(results) <- NULL

fwrite(results, file=out.fits.csv)

if ('snakemake' %in% ls()) {
    sink()
}
