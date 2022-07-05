#!/bin/bash

outdir=/n/data1/hms/dbmi/park/jluquette/glia/analysis/cov_heatmaps

muttypes=$(for x in $(ls -d try3/enrichment/boca/quantile/*___*); do basename $x; done)


for mut in $muttypes; do
    for nq in 3 5 10 50; do
        fs=$(ls try3/enrichment/*/quantile/${mut}/${nq}quantiles.csv)
        #echo "$mut $nq"
        #echo scripts/compute_enrichment_slopes.R \
            #$outdir/${mut}___${nq}quantiles___slopes.csv $fs
    done
done

# signature enrichment only exists for passA calls and for 3,5,10 quantiles
# indels currently not done?
for mut in neuron___A oligo___A; do #neuron___indel_A oligo___indel_A; do
    for nq in 3 5 10; do
        fs=$(ls try3/enrichment/*/quantile/${mut}/${nq}quantiles_sigenrich_adapted.csv)
        echo "$mut $nq"
        scripts/compute_enrichment_slopes.R \
            $outdir/${mut}___${nq}quantiles___slopes___sigenrich.csv $fs
    done
done
