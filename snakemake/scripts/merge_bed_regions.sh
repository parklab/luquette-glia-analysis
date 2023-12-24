#!/bin/bash

if [ $# -ne 3 ]; then
    echo "usage: input.bed hg19.genome output.bed"
    echo "hg19.genome is used to order chromosomes"
    exit 1
fi

inbed=$1
genome=$2
outbed=$3

if [ -f $outbed ]; then
    echo "output file $outbed already exists, please delete it first"
    exit 1
fi


(grep '^#' $inbed ; \
 for region in $(cut -f4 $inbed |uniq|sort|uniq|grep -v '^#'); do \
        echo "region=$region" 1>&2
        echo "lines before: $(awk '$0 ~ /^#/ || $4 == "'$region'"' $inbed|wc -l)" 1>&2
        echo "lines after: $(awk '$0 ~ /^#/ || $4 == "'$region'"' $inbed | bedtools merge | awk 'BEGIN { OFS="\t"; } { print $0, "'$region'" }'|wc -l)" 1>&2
        awk '$0 ~ /^#/ || $4 == "'$region'"' $inbed | bedtools merge | awk 'BEGIN { OFS="\t"; } { print $0, "'$region'" }'; done | bedtools sort -g $genome) > $outbed
