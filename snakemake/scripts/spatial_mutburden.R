#!/usr/bin/env Rscript

# detect script being run by snakemake
# if so, make a mock commandArgs function
if ('snakemake' %in% ls()) {
    logfile <- snakemake@log[[1]]
    con <- file(logfile, 'w')
    sink(con, type='output')
    sink(con, type='message')

    commandArgs <- function(...) unlist(c(
        snakemake@input['muts'],
        snakemake@input['meta'],
        snakemake@input['binbed'],
        snakemake@params['rolling_window_size'],
        snakemake@output['profile'],
        snakemake@output['stats'],
        snakemake@output['pdf'],
        snakemake@output['svg']
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 8) {
    #cat('typically power is lower for finding regions with significant depletion.\n')
    #cat('excluding these regions from consideration a priori allows fewer hypothesis tests.\n')
    stop('usage: spatial_mutburden.R muts_FILTERED.csv metadata.csv tiles.bed rolling_window_size out.profile.csv out.stats.csv out.modelcomparison.pdf out.modelcomparison.svg')
}


in.muts <- args[1]
in.meta.csv <- args[2]
in.bin.bed <- args[3]
bins.per.window <- as.integer(args[4])
out.profile.csv <- args[5]
out.stats.csv <- args[6]
out.pdf <- args[7]
out.svg <- args[8]

for (f in c(out.profile.csv, out.stats.csv)) {
    if (file.exists(f))
        stop(paste('output file', f, 'already exists, please delete it first'))
}

suppressMessages(library(scan2))
suppressMessages(library(mutenrich))
suppressMessages(library(MASS))
suppressMessages(library(extrafont))
suppressMessages(library(svglite))
if (!("Arial" %in% fonts()))
    stop("Arial font not detected; did you load extrafonts and run font_import() with the appropriate path?")

padj.method <- 'BH'

muts <- fread(in.muts)
gmuts <- mutenrich::gr(muts)

# To get sum(age) for mutations per mb*year
meta <- fread(in.meta.csv)

# In years
missing.samples <- setdiff(unique(muts$sample), meta$sample)
if (length(missing.samples) > 0) {
    warning(paste('metadata file', in.meta.csv, 'does not contain all samples in mutation table', in.muts, ':', paste(missing.samples, collapse=', ')))
}
agesum <- sum(meta[sample %in% unique(muts$sample)]$age)

binbed <- fread(in.bin.bed)
colnames(binbed) <- c('chr', 'start', 'end', 'id', 'keep', 'mean.depth')
binbed <- binbed[keep == 1]
bins <- GRanges(seqnames=sub('chr', '', binbed$chr), ranges=IRanges(start=binbed$start, end=binbed$end))

bins$muts <- countOverlaps(bins, gmuts)

# stats::filter - although oddly named, this computes a moving average
# using the `filter` argument as coefficients. coef=1 -> moving sum
# when the length of `filter` is even, the extra window is on the
# right.
# Result: sum of mutations over a sliding window of size 10 bins.
bins$muts.smoothed <- as.numeric(stats::filter(bins$muts, filter=rep(1, bins.per.window), sides=2))
bins$muts.smoothed.peryear <- bins$muts.smoothed / agesum


# Use two models to fit a constant background rate: Poisson and
# negative binomial (to handle a slight overdispersion, which was
# found by manual data exploration.
#
# Using GLM to fit MLEs. Mostly to handle negative binomial fitting.
# Mean values are log(base e).
pois.model <- glm(muts.smoothed ~ 1, data=mcols(bins), family=poisson)
negbin.model <- MASS::glm.nb(muts.smoothed ~ 1, data=mcols(bins))

# IMPORTANT: these are not independent tests because a sliding window is used!
bins$pois.pval <- sapply(bins$muts.smoothed, function(n) {
    if (is.na(n)) NA else poisson.test(n, r=exp(coef(pois.model)[1]))$p.value
})
bins$pois.pval.adj <- p.adjust(bins$pois.pval, method=padj.method)


bins$pois.pval.enrich <- sapply(bins$muts.smoothed, function(n) {
    mu <- exp(coef(pois.model)[1])
    if (is.na(n) | n <= mu) NA else poisson.test(n, r=mu, alternative='greater')$p.value
})
bins$pois.pval.enrich.adj <- p.adjust(bins$pois.pval.enrich, method=padj.method)

# Two-sided negative binomial test. This fragment was copied from poisson.test
# and all instances of {d,p}pois replaced by {d,p}nbinom.  Still need to
# carefully review to make sure it makes sense.
negbin.test <- function(x, mu, theta, relErr=1 + 1e-7) {
    if (is.na(x))
        return(NA)

    d <- dnbinom(x, mu=mu, size=theta)
    if (x == mu) {
        1 
    } else if (x < mu) {
        N <- ceiling(2 * mu - x)
        while (dnbinom(N, mu=mu, size=theta) > d)
            N <- 2 * N
        i <- seq.int(from = ceiling(mu), to = N)
        y <- sum(dnbinom(i, mu=mu, size=theta) <= d * relErr)
        pnbinom(x, mu=mu, size=theta) + pnbinom(N - y, mu=mu, size=theta, lower.tail = FALSE)
    } else {
        i <- seq.int(from = 0, to = floor(mu))
        y <- sum(dnbinom(i, mu=mu, size=theta) <= d * relErr)
        pnbinom(y - 1, mu=mu, size=theta) + pnbinom(x - 1, mu=mu, size=theta, lower.tail = FALSE)
    }
}


bins$negbin.pval <- sapply(bins$muts.smoothed, function(n) {
    negbin.test(n, mu=exp(coef(negbin.model)[1]), theta=negbin.model$theta)
})
bins$negbin.pval.adj <- p.adjust(bins$negbin.pval, method=padj.method)


negbin.test.enrich <- function(x, mu, theta)
    pnbinom(x - 1, mu=mu, size=theta, lower.tail=FALSE)

# Only test windows with > mean signal.
bins$negbin.pval.enrich <- sapply(bins$muts.smoothed, function(n) {
    mu <- exp(coef(negbin.model)[1])
    if (!is.na(n) & n > mu)
        negbin.test.enrich(n, mu=mu, theta=negbin.model$theta)
    else NA
})
bins$negbin.pval.enrich.adj <- p.adjust(bins$negbin.pval.enrich, method=padj.method)

# Write some select summary statistics to a CSV file.
writeLines(text=c(paste0('AgeSum,', agesum),
  paste0('Missing samples for agesum,', paste(missing.samples, collapse=',')),
  paste0('N cells,', length(unique(muts$sample))),
  paste0('N bins(alignable),', length(bins)),
  paste0('Bin width(bp),', width(bins)[1]),
  paste0('Bins per sliding window,', bins.per.window),
  paste0('N mutations,', nrow(muts)),
  paste0('Mutation type included,', unique(muts$muttype)),
  paste0('P-value adjustment method,', padj.method),
  paste0('Poisson GLM mean,', exp(coef(pois.model))),
  paste0('NegBin GLM mean,', exp(coef(negbin.model))),
  paste0('NegBin GLM theta(dispersion),', negbin.model$theta)),
    con=out.stats.csv)

out.df <- as.data.frame(bins)
colnames(out.df)[1] <- 'chr'  # replace 'seqnames'
out.df$strand <- NULL
setDT(out.df)
fwrite(out.df, file=out.profile.csv)

devs=list(pdf, svglite)
outs=c(out.pdf, out.svg)
for (i in 1:2) {
    devs[[i]](width=4.0, height=1.5, pointsize=5, file=outs[i])
    par(mar=c(4,4,1,0))
    layout(matrix(c(1:3,4,4,4), nrow=2, byrow=TRUE), widths=c(2,1,2))
    xmax <- max(bins$muts.smoothed, na.rm=TRUE)
    plot(as.numeric(names(table(bins$muts.smoothed))),
         as.vector(table(bins$muts.smoothed)/length(bins)),
        type='h', log='y', bty='n',
        ylab='Fraction of windows', xlab='Mutations per window')
    # plot extra points on top of histogram bars
    points(as.numeric(names(table(bins$muts.smoothed))),
         as.vector(table(bins$muts.smoothed)/length(bins)), pch=20)
    lines(0:xmax, dpois(0:xmax, lambda=exp(coef(pois.model)[1])), lwd=1, col=1)
    lines(0:xmax, dnbinom(0:xmax, mu=exp(coef(negbin.model)[1]), size=negbin.model$theta), col=2, lwd=1)
    par(mar=c(0,0,0,0))
    plot(0, pch=NA, xlab='', ylab='', bty='n', xaxt='n', yaxt='n')
    legend('topleft', col=1:2, lwd=2, legend=c('Poisson', 'Neg. Bin.'), bty='n', title='Fit method')
    par(mar=c(4,4,1,1))
    plot(-log10(bins$negbin.pval), -log10(bins$pois.pval), bty='n',
        ylab='-log10(Poisson p-value)', xlab='-log10(NegBin p-value)')
    abline(coef=0:1)

    # TODO: generalize these power calculations so they're done for every analysis
    # XXX: add enrich-only power testing
    if (FALSE) {
    mu <- exp(coef(pois.model))
    mults=c(1,2,5,10,20,50,100)

    for (i in 1:length(mults)) {
        mult = mults[i]
        curve(ppois(q=qpois(p=0.025/length(bins), lambda=mu*mult), lambda=x*mu*mult),
            from=0.1,to=1,log='x', main='Power vs. effect size for depletion\n1MB smoothed at 100kb, Bonferroni', ylab='Power', xlab='Effect size', add=i>1, col=i, lwd=2, ylim=c(0,1))
    }

    for (i in 1:length(mults)) {
        mult = mults[i]
        curve(1-ppois(q=qpois(p=1-0.025/length(bins), lambda=mu*mult, lower.tail=T), lambda=x*mu*mult), from=1,to=10,log='x', main='Power vs. effect size for enrichment\n1MB smoothed at 100kb, Bonferroni', ylab='Power', xlab='Effect size', add=i>1, col=i, lwd=2, ylim=c(0,1))
    }
    legend('bottomright', legend=mults, col=1:length(mults),lwd=2, title='Mutation multiplier')
    dev.off()
    }
}

if ('snakemake' %in% ls()) {
    sink()
}

