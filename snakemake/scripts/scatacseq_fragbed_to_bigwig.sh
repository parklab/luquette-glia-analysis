#!/bin/bash

set -euo pipefail

if [ $# -ne 4 ]; then
    echo "usage: $0 fragments.tsv hg19.genome out.bigwig tmpfile"
    exit 1
fi

infragbed=$1
genome=$2
outbigwig=$3
tmpfile=$4

if [ -f $outbigwig ]; then
    echo "output file $outbigwig already exists, please delete it first"
    exit 1
fi
if [ -f $tmpfile ]; then
    echo "output file $tmpfile already exists, please delete it first"
    exit 1
fi


# bedtools sort doesn't work for some reason. bedGraphToBigWig
# will complain with:
# /dev/stdin is not case-sensitive sorted at line 754830.  Please use "sort -k1,1 -k2,2n" with LC_COLLATE=C,  or bedSort and try again.
# So I guess we just sort as it suggests.
#| bedtools sort -i /dev/stdin -g $genome ` \
echo "running genomecov"
bedtools genomecov -i $infragbed -g $genome -bga \
    | sort -k1,1 -k2,2n > $tmpfile 

echo "running bedGraphToBigWig"
bedGraphToBigWig $tmpfile $genome $outbigwig
