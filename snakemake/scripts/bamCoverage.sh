#!/bin/bash
#SBATCH -p short
#SBATCH -A park
#SBATCH -t 12:00:00
#SBATCH --mem=4G

if [ $# -ne 2 ]; then
    echo "usage: $0 in.bam out.bw"
    exit 1
fi

inbam=$1
outbw=$2

if [ -f $outbw ]; then
    echo "output file $outbw already exists, please delete it first"
    exit 1
fi

bamCoverage -v \
    --binSize 100 \
    --minMappingQuality 60 \
    --samFlagExclude 3840 \
    -b $inbam \
    -o $outbw
