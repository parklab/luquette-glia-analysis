#!/bin/bash

set -eo pipefail

if [ $# -ne 2 ]; then
    echo "usage: $0 input.vcf output.vcf"
    exit 1
fi

snpeff=snakemake/scripts/snpEff/snpEff.jar
config=snakemake/scripts/snpEff/snpEff.config

invcf=$1
outvcf=$2

if [ -f $outvcf ]; then
    echo "output file $outvcf already exists, please delete it first"
    exit 1
fi


java -jar $snpeff -c $config -t -noStats -v hg19 $invcf > $outvcf
