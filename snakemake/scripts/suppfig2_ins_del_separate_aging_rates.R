#!/usr/bin/env Rscript

# detect script being run by snakemake
# if so, make a mock commandArgs function
if ('snakemake' %in% ls()) {
    logfile <- snakemake@log[[1]]
    con <- file(logfile, 'w')
    sink(con, type='output')
    sink(con, type='message')

    commandArgs <- function(...) unlist(c(
        snakemake@input['meta'],
        snakemake@input['muts'],
        snakemake@output['data_csv'],
        snakemake@output['pdf'],
        snakemake@output['models_csv'],
        snakemake@params['color']
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}


args <- commandArgs(trailingOnly=T)

if (length(args) != 6)
    stop(sprintf("usage: suppfig_ins_del_separate_aging_rates.R sample_metadata.csv muts.csv out.data.csv out.pdf out.model.csv plot_color"))

meta.file <- args[1]
muts.file <- args[2]
out.data.csv <- args[3]
out.pdf <- args[4]
out.models.csv <- args[5]
color <- args[6]

for (f in c(out.data.csv, out.models.csv, out.pdf)) {
    if (file.exists(f))
        stop(sprintf("output file %s already exists, please delete it first", f))
}

# DO WORK
suppressMessages(library(data.table))
suppressMessages(library(lme4))
suppressMessages(library(lmerTest))

meta <- fread(meta.file)[,.(donor, sample, age)]

# don't filter down to indels first, samples with 0 indels will drop out
muts <- fread(muts.file)#[muttype == 'indel']
muts[, indel.type := ifelse(muttype == 'indel', ifelse(nchar(altnt) - nchar(refnt) > 0, 'ins', 'del'), 'not_indel')]

id.tab <- meta[muts[,.(ins=sum(indel.type=='ins'), del=sum(indel.type=='del')),by=.(sample)],,on=.(sample)]
fwrite(id.tab, file=out.data.csv)

ins.model <- lmer(ins ~ age + (1|donor), data=id.tab)
del.model <- lmer(del ~ age + (1|donor), data=id.tab)

sum.to.dt <- function(m, indel.type)
    cbind(indel.type=indel.type, var=c('Intercept', 'Age'), as.data.table(coef(summary(m))))

fwrite(rbind(sum.to.dt(ins.model, 'Insertion'), sum.to.dt(del.model, 'Deletion')), file=out.models.csv)


pdf(file=out.pdf, width=2.5, height=2.5, pointsize=5)
plot(id.tab[,.(age, del)], type='p', pch=17, col=color, xlab='Age', ylab='Somatic indel calls')
points(id.tab[,.(age, ins)], type='p', pch=1, cex=1, lwd=1.2, col=color)

# don't plot below 0
m <- ins.model
curve(coef(summary(m))[1] + coef(summary(m))[2]*x,
            from=0, to=2*max(id.tab$age), lwd=1, col=color, add=T, lty='dashed')
m <- del.model
curve(coef(summary(m))[1] + coef(summary(m))[2]*x,
            from=0, to=2*max(id.tab$age), lwd=1, col=color, add=T, lty='solid')

legend("topleft", pch=c(17,1), lty=c('solid', 'dashed'), legend=c('Deletion', 'Insertion'), bty='n', col=color)

dev.off()


if ('snakemake' %in% ls()) {
    sink()
}
