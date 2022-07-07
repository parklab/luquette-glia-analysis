#!/usr/bin/env python
#SBATCH -p park
#SBATCH -A park_contrib
#SBATCH -t 24:00:00
#SBATCH --mem=4G

from argparse import ArgumentParser

ap = ArgumentParser()
ap.add_argument('fragments_tsv',
    help='fragments.tsv file as output by 10x CellRanger software.')
ap.add_argument('sample_to_cell_map',
    help='Tab delimited text file with two columns: column 1=cell tag and library ID (e.g., AAACGAAAGAAATACC-9); column 2=cell type (e.g., excitatory neuron).')
ap.add_argument('output_prefix',
    help='Prefix to add to output files. There will be one output file for each (library ID, cell type) combination and one merged output file for each cell type.')
ap.add_argument('-n', metavar='INT', default=11,
    help='Number of libraries in merged output. 11 libraries were used in this study.')
args = ap.parse_args()

tag_to_celltype = {}
with open(args.sample_to_cell_map, 'r') as infile:
    for line in infile:
        tag, celltype = line.strip().split('\t')
        celltype = celltype.replace(' ', '_')
        tag_to_celltype[tag] = celltype

print('got ' + str(len(tag_to_celltype)) + ' cells spanning ' + str(len(set(tag_to_celltype.values()))) + ' cell types')

mergedfiles = [ (ct, open(args.output_prefix + '.librarymerged_' + ct + '.bed', 'w')) for ct in set(tag_to_celltype.values()) ]
perlibfiles = [ (str(i) + ' ' + ct, open(args.output_prefix + '.library' + str(i) + '_' + ct + '.bed', 'w'))
                 for i in range(1, args.n+1) for ct in set(tag_to_celltype.values()) ]

outfiles = dict(mergedfiles + perlibfiles)
#print(outfiles)

with open(args.fragments_tsv, 'r') as infile:
    lineno=0
    failedlines=0
    for line in infile:
        lineno = lineno+1
        chrom, start, stop, tag, _ = line.strip().split('\t')
        start = int(start)
        stop = int(stop)

        # Both the start and stop of the fragment represent transposition
        # events, so both need to be output as separate signals.

        # Format is: <barcode>-<library number>, e.g. GCAGATTCAGTGCGTC-6
        lib = tag.split('-')[1]
        celltype = tag_to_celltype.get(tag, None)
        if celltype is not None:
            of = outfiles[str(lib) + ' ' + celltype]
            of.write(chrom + '\t' + str(start-1) + '\t' + str(start) + '\n')
            of.write(chrom + '\t' + str(stop-1) + '\t' + str(stop) + '\n')
            of = outfiles[celltype]
            of.write(chrom + '\t' + str(start-1) + '\t' + str(start) + '\n')
            of.write(chrom + '\t' + str(stop-1) + '\t' + str(stop) + '\n')
        else:
            failedlines = failedlines+1


print('lines:         ' + str(lineno))
print('failed lines:  ' + str(failedlines))
for outf in outfiles.values():
    print('closing file: ' + str(outf))
    outf.close()
