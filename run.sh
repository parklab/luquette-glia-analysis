#!/bin/bash

dryrun=$1

if [ "x$dryrun" != 'x' ]; then
    flags="--dryrun --quiet"
fi

snakemake $flags \
    --dir try3 \
    -s Snakefile \
    -C INPUT_MANIFEST=`realpath INPUT_MANIFEST` ROADMAP_EIDS=`realpath ROADMAP_EIDS_WITH_H3K27ac`
