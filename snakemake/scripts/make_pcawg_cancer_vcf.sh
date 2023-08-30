#!/bin/bash

set -eo pipefail

if [ $# -lt 2 ]; then
    echo "usage: $0 tumor_type vcf1 [ vcf2 ... vcfN ]"
    exit 1
fi


tumor=$1
shift 1
files=$@

echo '##fileformat=VCFv4.1'
echo '##INFO=<ID=TType,Number=1,Type=String,Description="Tumor or cell type">'
echo '##INFO=<ID=Origin,Number=1,Type=String,Description="File of origin">'
echo '##INFO=<ID=Sample,Number=1,Type=String,Description="Sample of origin">'
for f in $files; do
    bcftools view -v snps $f 2> /dev/null \
        | vcfkeepinfo - . \
        | grep -v '#' \
        | awk 'BEGIN { OFS="\t"; }{ print $1,$2,$3,$4,$5,$6,$7,"TType='$tumor';Origin='$(basename $f)';Sample=."}'
done | sort -V
