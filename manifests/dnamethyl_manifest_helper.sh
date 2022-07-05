#!/bin/bash

echo "datasource,dataclass,experiment_type,signal_type,eid,file"
for x in data/roadmap/dnamethyl/??BS_*; do
    a=$(basename $x|cut -f1 -d_)
    b=$(basename $x|cut -f2 -d_)
    for f in $x/bigwig/*.bigwig; do
        eid=$(basename $f|cut -c1-4)
        echo encode,dnamethyl,$a,$b,$eid,$f
    done
done
