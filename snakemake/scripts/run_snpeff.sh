#!/bin/bash

if [ $# -ne 2 ]; then
    echo "usage: $0 input.vcf output.vcf"
    exit 1
fi

snpeff=snakemake/scripts/snpEff.jar

invcf=$1
outvcf=$2

if [ -f $outvcf ]; then
    echo "output file $outvcf already exists, please delete it first"
    exit 1
fi


java -jar $snpeff -t -noStats -v hg19 $invcf > $outvcf
