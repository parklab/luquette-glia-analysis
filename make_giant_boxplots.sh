#!/bin/bash

outdir=/n/data1/hms/dbmi/park/jluquette/glia/analysis/giant_boxplots

for f in try3/enrichment/*/quantile/*/*quantiles.csv; do
    outbase=$(echo $f | cut -f 3- -d/ | tr '/' '_' | sed -e's/.csv$//')
    echo "$f -> $outbase"
    scripts/boxplot_enrichment_grid.R $f \
            enrichment_grid_R_ignore_none \
            $outdir/$outbase.svg \
            $outdir/$outbase.pdf \
            $outdir/$outbase.csv
done
