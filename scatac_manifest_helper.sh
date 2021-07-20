#!/bin/bash

# manually encoded library<N> -> sample name map via sed
for f in ../figures/scatacseq/bigwig/*.bigwig; do
    libstr=$(echo $f|cut -f4 -d\.)
    libid=$(echo $libstr | cut -f1 -d_)
    celltype=$(echo $libstr|cut -f2- -d_)
    echo scatacseq,$libid,$celltype,data/$(echo $f|cut -f3- -d/)
done | sed -e 's/,library1,/,library1,1278,/' -e's/,library2,/,library2,1465-S3PFC-1_03132020_Analysis,/' -e 's/,library3,/,library3,1465-S3PFC-2_03132020_Analysis,/' -e's/,library4,/,library4,4638,/' -e 's/,library5,/,library5,4643,/' -e's/,library6,/,library6,5087_highdp,/' -e's/,library7,/,library7,5219_highdp,/' -e 's/,libr
