#!/bin/bash

set -eo pipefail

if [ $# -ne 4 ]; then
    echo "usage: donor_id input_maf output_maf output_dir"
    echo "you must have vcf2maf.pl in your path (conda activate scan2)"
    exit 1
fi

did=$1
inmaf=$2
outmaf=$3
outvcf=${outmaf/.maf/.vcf}
outdir=$4

if [ -f $outmaf ]; then
    echo "output MAF file $outmaf already exists, please delete it first"
    exit 1
fi

if [ "x$(echo $outmaf | grep '.maf$')" != "x$outmaf" ]; then
    echo "$outmaf must end in .maf"
    exit 1
fi

if [ -f $outvcf ]; then
    echo "output VCF file $outvcf (derived from output_maf) already exists, please delete it first"
    exit 1
fi

(head -1 $inmaf ; awk  -F$'\t' '$43 == "'$did'"' $inmaf) > $outmaf

maf2vcf.pl --input-maf $outmaf \
    --output-dir $outdir \
    --output-vcf $outvcf \
    --ref-fasta ~/ndata1/genotyper1/paper/resources/human_g1k_v37_decoy.fasta \
