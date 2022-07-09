#!/bin/bash

echo 'datasource,celltype,geoid,file'
wget --quiet -O /dev/stdout http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeUwRepliSeq/files.txt|grep WaveSignal | while read line; do f=$(echo $line|cut -f1 -d\ ); cell=$(echo $line|cut -f7 -d\;|sed -e 's/ cell=//'); geo=$(echo $line|grep -o 'geoSampleAccession=[^;]*'|cut -f2 -d=); echo repliseq,$cell,$geo,data/repliseq/bigwig/$f; done|sort -k3
