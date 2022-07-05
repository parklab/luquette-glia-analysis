#!/bin/bash

echo "datasource,dataclass,model,eid,file"

# 15state classifications
for f in ../../pta/chromhmm/beds_all_tissues15/E*_mnemonics.bed; do
    eid=$(basename $f|cut -c1-4)
    echo "roadmap,chromhmm,15state,$eid,data/roadmap/chromhmm/15state/bed/${eid}___200bp_tiles.bed"
done

# 18state classifications
for f in ../../pta/chromhmm/beds_all_tissues18/E*_mnemonics.bed; do
    eid=$(basename $f|cut -c1-4)
    echo "roadmap,chromhmm,18state,$eid,data/roadmap/chromhmm/18state/bed/${eid}___200bp_tiles.bed"
done
