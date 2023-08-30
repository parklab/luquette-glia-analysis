#!/bin/bash
#SBATCH -p park
#SBATCH -A park_contrib
#SBATCH -t 12:00:00
#SBATCH --mem=7000

set -eo pipefail

if [ $# -ne 3 ]; then
    echo "usage: $0 in.bigwig temporar.out.bedgraph out.bigwig"
    exit 1
fi

inbw=$1
outbg=$2
outbw=$3

if [ -f $outbg ]; then
    echo "output file $outbg already exists, please delete it first"
    exit 1
fi

if [ -f $outbw ]; then
    echo "output file $outbw already exists, please delete it first"
    exit 1
fi

echo "Making bedgraph $outbg.."
bigWigToBedGraph $inbw /dev/stdout \
    | grep -v '^GL' | grep -v '^NC' | grep -v '^hs' \
    | sed -e's/^/chr/' -e 's/MT/M/' \
    | sort -k1,1 -k2,2n > $outbg

echo "Converting bedgraph to final bigwig $outbw.."
bedGraphToBigWig $outbg /n/data1/hms/dbmi/park/jluquette/pta/spatial_depth/hg19.chrom.sizes $outbw
