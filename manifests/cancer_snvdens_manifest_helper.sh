#!/bin/bash

echo "datasource,tumor,signal_type,file"
cat /n/data1/hms/dbmi/park/jluquette/glia/icgc/tumor_types.txt | while read tumor; do
    echo cancer_snvdens,$tumor,sumdens,data/cancer_snvdens/bigwig/${tumor}_sumdens.bigwig
    echo cancer_snvdens,$tumor,normdens,data/cancer_snvdens/bigwig/${tumor}_normdens.bigwig
done
