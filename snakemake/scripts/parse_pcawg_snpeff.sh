#!/bin/bash

if [ $# -ne 1 ]; then
    echo "usage: $0 input.vcf"
    exit 1
fi

invcf=$1

cat $invcf \
    | grep -v '#' \
    | grep ANN \
    | cut -f 8 \
    | tr ';' '\t' \
    | sed 's#TType=##g' \
    | sed 's#Origin=##g' \
    | sed 's#Sample=##g' \
    | sed 's#ANN=##g' \
    | cut -f 1,4 \
    | awk 'BEGIN { OFS="\t"; } { n = split($2, anns, ","); for (i = 1; i <= n; ++i) { print NR, $1, anns[i];  } }' \
    | tr '|' '\t' \
    | awk 'BEGIN { OFS = "\t"; } { print $1, $2, $6 "/" $7; }'
