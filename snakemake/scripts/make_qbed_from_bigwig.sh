#!/bin/bash

set -euo pipefail

if [ $# -lt 6 ]; then
    echo "usage: $0 n_quantiles min_coverage tiles.bed signal.{bigwig|sig} out.qbed tmpout tmpout2 [ metadata1 ... metadataN ]"
    echo "min_coverage - a fraction in the range [0,1]. Tiles for which the bigwig"
    echo "    signal covers < min_coverage bases (as a fraction of tile size) will"
    echo "    be assigned an NA (excluded) score and quantile."
    echo "    Set to 0 to disable this functionality and behave like the old script."
    echo "if the signal file ends in the (case sensitive) .sig extension, then"
    echo "    bigWigAverageOverBed is not run it is assumed that the .sig file"
    echo "    is in the correct format."
    echo "metadata should be supplied as TAG=VALUE pairs. Neither TAG nor VALUE"
    echo 'may contain a semi-colon (;).'
    exit 1
fi

nquantiles=$1
mincov=$2
tilesbed=$3
bigwig="$4"
outqbed="$5"
tmpout="$6"
tmpout2="$7"
shift 7
metadata=("QUANTILES=$nquantiles" "MIN_COVERAGE=$mincov" "$@")

echo "Minimum coverage requirement: $mincov"

if [ "x$(echo "$bigwig" | grep -c '\.sig$')" == "x0" ]; then
    is_sig=FALSE
    echo "Signal file is type=BIGWIG"
else
    is_sig=TRUE
    echo "Signal file is type=SIGNAL, a file produced by bigWigAverageOverBed"
fi
    
if [ "x$(echo "$metadata" | grep -c ';')" != "x0" ]; then
    echo "metadata TAG=VALUE pairs must not contain semicolons"
    exit
fi

metatags="#QBED_VERSION=1"
for kv in "${metadata[@]}"; do
    if [ "x$(echo "$kv"|grep -c '=')" == "x0" ]; then
        echo "key-value pair '$kv' does not contain an = sign"
        exit 1
    fi
    metatags="$metatags;$kv"
done
echo "using metatags=$metatags"


if [ -f "$outqbed" ]; then
    echo "output file $outqbed already exists, please delete it first"
    exit 1
fi
if [ -f "$tmpout" ]; then
    echo "output file $tmpout already exists, please delete it first"
    exit 1
fi
if [ -f "$tmpout2" ]; then
    echo "output file $tmpout2 already exists, please delete it first"
    exit 1
fi


# column 5 of tmpout is mean with non-covered bases counted as 0
# column 6 of tmpout is mean over just covered bases
if [ $is_sig == "FALSE" ]; then
    echo "Averaging bigwig over tiles.."
    ( bigWigAverageOverBed "$bigwig" $tilesbed "$tmpout" )
else
    tmpout="$bigwig"
    echo "Skipping bigWigAverageOverBed because type=SIGNAL, using tmpout=$tmpout"
fi


# Map scores to quantiles and write out the quantiles for every tile,
# including NA quantiles for excluded tiles.
# Column 5 of the tile bed is 1 if the tile is considered analyzable
# and 0 otherwise.
echo "Converting to quantiles (n=$nquantiles)"
Rscript -e 'library(data.table); keep <- fread("'$tilesbed'")[[5]]; tmpout <- fread("'$tmpout'"); score <- tmpout[[5]]; coverage <- tmpout[[3]]/tmpout[[2]]; score[keep == 0 | coverage < '$mincov'] <- NA; q <- findInterval(score, quantile(score, na.rm=T, probs=1:'$nquantiles'/'$nquantiles'), rightmost.closed=TRUE)+1; cat("gc before writing\n"); print(gc()); write(paste(q, score, sep="\t"), ncolumns=1, file="'$tmpout2'"); cat("gc after writing\n"); print(gc())'

# tilesbed has 5 columns, tmpout2 has 2
echo "Writing qbed '$outqbed'.."
(echo "$metatags" ;
 paste $tilesbed "$tmpout2" | cut -f1-3,6-7) > "$outqbed"
