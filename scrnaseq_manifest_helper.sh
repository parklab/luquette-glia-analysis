#!/bin/bash

echo "datasource,dataclass,donor,selection,celltype,file"
cut -f2- -d, ../scrnaseq/rds/combined_mean_expr.csv |head -1 | tr ',' '\n' | while read line; do
    donor=$(echo $line | cut -f1 -d\.)
    suffix=$(echo $line | cut -f2 -d\. | sed -e 's/___/ /')
    selection=$(echo $suffix | cut -f1 -d\ )
    celltype=$(echo $suffix | cut -f2 -d\ )
    echo scrnaseq,expression,$donor,$selection,$celltype,data/scrnaseq/expression/bigwig/$line.bigwig
done

# Manual output for the cross-library combined celltypes
for ct in Astrocytes Microglia Endothelial Oligodendrocytes OPCs Neurons Inhibitory-Neurons Excitatory-Neurons; do
    donor=combined
    selection=combined
    celltype=$ct
    echo scrnaseq,expression,$donor,$selection,$celltype,data/scrnaseq/expression/bigwig/${ct}.bigwig
done
